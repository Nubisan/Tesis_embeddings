"""
CSCLP (Constrained-Size Clustering via Linear Programming).
"""

from __future__ import annotations

import logging
import time
from typing import Any, Optional, cast

import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from scipy.sparse import vstack as sparse_vstack
from scipy.spatial.distance import pdist, squareform

from ._common import (
    compute_metrics,
    get_predictions_dir,
    pam_kaufman_rousseeuw,
    prepare_data,
    tabulate_clusters,
    to_int_list,
)

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(message)s")


CSCLP_HYPERPARAMS_GLOBAL: dict[str, Any] = {
    "algorithm":       "CSCLP",
    "distance_method": "cosine",
    "seed":            123,
    "lp_solver":       "scipy.optimize.linprog (HiGHS dual simplex)",
    "all_bin":         True,
}

CSCLP_RUNS: dict[str, dict[str, Any]] = {}


def _build_ilp_constraints(
    k: int, m: int, target: np.ndarray
) -> tuple[csr_matrix, np.ndarray]:
    """
    Restricciones del modelo en formato sparse.

    Variables: x[j*m + i] = 1 si el no-medoid i va al cluster j (j en [0,k), i en [0,m)).
    R1 (m filas): cada punto va a exactamente 1 cluster  -> sum_j x[j*m+i] = 1
    R2 (k filas): cada cluster j recibe target[j]-1 puntos -> sum_i x[j*m+i] = target[j]-1
    """
    n_vars = k * m

    # R1: cada punto en un solo cluster
    rows_p = np.tile(np.arange(m), k)
    cols_p = np.repeat(np.arange(0, k * m, m), m) + np.tile(np.arange(m), k)
    data_p = np.ones(m * k, dtype=float)
    a_point = csr_matrix((data_p, (rows_p, cols_p)), shape=(m, n_vars))
    b_point = np.ones(m, dtype=float)

    # R2: cardinalidad por cluster
    rows_c = np.repeat(np.arange(k), m)
    cols_c = np.arange(k * m)
    data_c = np.ones(k * m, dtype=float)
    a_card = csr_matrix((data_c, (rows_c, cols_c)), shape=(k, n_vars))
    b_card = (target - 1).astype(float)

    if np.any(b_card < 0):
        raise ValueError("target_cardinality[j] - 1 < 0, modelo inválido.")

    a_stacked_any: Any = sparse_vstack([a_point, a_card])
    a_full = cast(csr_matrix, a_stacked_any.tocsr())
    b_full = np.concatenate([b_point, b_card])
    return a_full, b_full


def run_csclp_algo(
    dataset: pd.DataFrame,
    target_cardinality: list[int],
    dist_method: str = "cosine",
    seed: int = 123,
) -> dict[str, Any]:
    """
    Ejecuta CSCLP: prepara X/y, calcula distancias N×N, PAM para k medoids y
    asignación de no-medoids por relajación lineal (óptimo entero por unimodularidad).
    """
    prepared = prepare_data(dataset)
    if prepared is None:
        raise ValueError("Dataset inválido (no se pudo preparar X/y).")

    x_df:     pd.DataFrame    = prepared["X"]
    y_factor: pd.Categorical = prepared["y"]

    x_mat  = x_df.to_numpy(dtype=float)
    n      = int(x_mat.shape[0])
    target = np.asarray(target_cardinality, dtype=int)
    k      = int(len(target))

    if k < 2:
        raise ValueError("k inválido (<2).")
    if k >= n:
        raise ValueError(f"k={k} inválido para n={n}.")
    if np.any(np.isnan(target.astype(float))):
        raise ValueError("target_cardinality contiene NaN.")
    if int(target.sum()) != n:
        raise ValueError(f"target_cardinality suma {int(target.sum())} pero n={n}.")
    if np.any(target < 1):
        raise ValueError("target_cardinality tiene valores < 1.")

    # El cronómetro abarca el pipeline completo (∆ -> PAM -> construcción -> LP),
    # igual que los otros 7 algoritmos; no solo el solve.
    start_total = time.perf_counter()

    # Distancia N×N (∆); se reutiliza para el silhouette.
    condensed = cast(np.ndarray, pdist(x_mat, metric=dist_method))  # pyright: ignore[reportArgumentType, reportCallIssue]
    full_dist = squareform(condensed)

    # Consistencia con set.seed(seed) del R; PAM (random.start=FALSE) no consume
    # RNG, pero mantenemos la llamada para no alterar el orden en cambios futuros.
    _ = np.random.default_rng(seed)
    medoid_idx = pam_kaufman_rousseeuw(full_dist, k)

    all_idx        = np.arange(n)
    non_medoid_idx = np.setdiff1d(all_idx, medoid_idx, assume_unique=True)
    m              = int(len(non_medoid_idx))

    # cost_mat[j, i] = distancia medoid j -> no-medoid i (k × m), aplanado row-major.
    cost_mat = full_dist[np.ix_(medoid_idx, non_medoid_idx)]
    f_obj = cost_mat.flatten("C")

    a_full, b_full = _build_ilp_constraints(k, m, target)

    # Matriz TU (cada variable en 1 restricción de asignación y 1 de cardinalidad):
    # la relajación continua ya da el óptimo entero; highs-ds devuelve un vértice.
    res = linprog(
        c      = f_obj,
        A_eq   = a_full,
        b_eq   = b_full,
        bounds = (0.0, 1.0),
        method = "highs-ds",
    )
    total_time = time.perf_counter() - start_total

    res_any = cast(Any, res)
    if not res_any.success or res_any.x is None:
        raise RuntimeError(
            f"LP infeasible/failed. status={res_any.status} message={res_any.message}"
        )

    # Entera por unimodularidad; redondeo por seguridad de floats.
    x_sol = np.round(np.asarray(res_any.x, dtype=float)).astype(int)
    assign_mat = x_sol.reshape(k, m)

    y_pred = np.zeros(n, dtype=int)
    y_pred[medoid_idx] = np.arange(1, k + 1)                    # medoids -> labels 1..k
    y_pred[non_medoid_idx] = np.argmax(assign_mat, axis=0) + 1  # argmax por columna

    return {
        "y_predict":      y_pred,
        "y":              y_factor,
        "x_mat":          x_mat,
        "full_dist":      full_dist,
        "execution_time": float(total_time),
        "hyperparams": {
            "seed":         seed,
            "dist_method":  dist_method,
            "k":            k,
            "n":            n,
            "lp_solver":    "scipy.optimize.linprog (HiGHS dual simplex)",
            "all_bin":      True,
        },
    }


def run_clustering_row(
    dataset: pd.DataFrame,
    target_cardinality: list[int],
    dataset_name: str,
) -> Optional[dict[str, Any]]:
    """Ejecuta CSCLP sobre un dataset y devuelve la fila de resultados."""
    result = run_csclp_algo(dataset, target_cardinality)

    y_pred:    np.ndarray     = result["y_predict"]
    y_factor:  pd.Categorical = result["y"]
    x_mat:     np.ndarray     = result["x_mat"]
    full_dist: np.ndarray     = result["full_dist"]
    total_t:   float          = result["execution_time"]

    y_int = (y_factor.codes + 1).astype(int)
    metrics = compute_metrics(y_int, y_pred, x_mat, dist_matrix=full_dist)

    CSCLP_RUNS[dataset_name] = result["hyperparams"]

    target           = np.asarray(target_cardinality, dtype=int)
    cardinality_pred = tabulate_clusters(y_pred, len(target))

    return {
        "name":               dataset_name,
        "n":                  len(y_pred),
        "k":                  int(len(np.unique(y_pred))),
        "y_predict":          " ".join(str(int(c)) for c in y_pred),
        "y_true":             " ".join(str(int(c)) for c in y_int),
        "target_cardinality": " ".join(str(int(c)) for c in target),
        "cardinality_pred":   " ".join(str(int(c)) for c in cardinality_pred),
        "Execution_Time":     total_t,
        "ARI":                metrics["ARI"],
        "AMI":                metrics["AMI"],
        "NMI":                metrics["NMI"],
        "Silhouette_mean":    metrics["Silhouette_mean"],
    }


def run(odatasets_unique: pd.DataFrame) -> None:
    """Loop principal"""
    results: list[dict[str, Any]] = []

    for i in range(len(odatasets_unique)):
        print(f"\n--- Running CSCLP for dataset at position: {i + 1} ---", flush=True)
        try:
            dataset      = odatasets_unique["dataset"].iloc[i]
            dataset_name = odatasets_unique["name"].iloc[i]
            target_raw   = odatasets_unique["class_distribution_vector"].iloc[i]

            target = to_int_list(target_raw)
            if dataset is None or target is None or len(target) < 2:
                print("Dataset o cardinalidad inválida. Skipping...")
                continue

            row = run_clustering_row(dataset, target, dataset_name)
            if row is not None:
                results.append(row)
                print(f"Time (total): {row['Execution_Time']:.4f} s")
            else:
                print("Dataset skipped.")

        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error processing dataset at position {i + 1}: {e}")

    if results:
        df_out = pd.DataFrame(results)
        out_path = get_predictions_dir() / "pred_CSCLP.csv"
        df_out.to_csv(out_path, index=False)
        print(f"\nCSCLP finalizado. Archivo guardado en: {out_path}")
    else:
        print("\nCSCLP finalizado sin resultados válidos.")