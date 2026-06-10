"""
algorithms/csclp.py
===================
CSCLP (Constrained-Size Clustering via Linear Programming).

Algoritmo en dos etapas:
  1. PAM (k-medoids) para encontrar k medoids representativos
  2. Programación Lineal Entera (ILP) para asignar los no-medoids al medoid
     más cercano respetando las restricciones de cardinalidad
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Optional, cast

import numpy as np
import pandas as pd
# scipy no provee stubs perfectos; silenciamos los warnings de tipos para
# imports que sabemos que funcionan en runtime.
from scipy.optimize import Bounds, LinearConstraint, milp  # type: ignore[reportMissingTypeStubs]
from scipy.sparse import csr_matrix  # type: ignore[reportMissingTypeStubs]
from scipy.sparse import vstack as sparse_vstack  # type: ignore[reportMissingTypeStubs]
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


# ============================================================================
# Hiperparámetros y rutas de salida
# ============================================================================

CSCLP_HYPERPARAMS_GLOBAL: dict[str, Any] = {
    "algorithm":       "CSCLP",
    "distance_method": "cosine",
    "seed":            123,
    "lp_solver":       "scipy.optimize.milp (HiGHS)",
    "all_bin":         True,
}

# Hiperparámetros guardados por dataset durante la ejecución
CSCLP_RUNS: dict[str, dict[str, Any]] = {}

# La ruta de salida se obtiene en runtime via get_predictions_dir()
# de _common, que lee la env var CLUSTERING_MODEL seteada por testing.py.


# ============================================================================
# Construcción del ILP
# ============================================================================

def _build_ilp_constraints(
    k: int, m: int, target: np.ndarray
) -> tuple[csr_matrix, np.ndarray]:
    """
    Construye las restricciones del ILP en formato sparse.

    Variables: x[j*m + i] = 1 si el no-medoid i va al cluster j.
        - j en [0, k), i en [0, m)
        - Total: n_vars = k * m

    Restricción 1 (m filas): cada punto i va a exactamente 1 cluster
        sum_j x[j*m + i] = 1     para todo i

    Restricción 2 (k filas): cada cluster j recibe target[j]-1 no-medoids
        sum_i x[j*m + i] = target[j] - 1   para todo j
    """
    n_vars = k * m

    # --- Restricción 1: A_point[i, j*m + i] = 1 ---
    rows_p = np.tile(np.arange(m), k)
    cols_p = np.repeat(np.arange(0, k * m, m), m) + np.tile(np.arange(m), k)
    data_p = np.ones(m * k, dtype=float)
    a_point = csr_matrix((data_p, (rows_p, cols_p)), shape=(m, n_vars))
    b_point = np.ones(m, dtype=float)

    # --- Restricción 2: A_card[j, j*m + i] = 1 para todo i ---
    rows_c = np.repeat(np.arange(k), m)
    cols_c = np.arange(k * m)
    data_c = np.ones(k * m, dtype=float)
    a_card = csr_matrix((data_c, (rows_c, cols_c)), shape=(k, n_vars))
    b_card = (target - 1).astype(float)

    if np.any(b_card < 0):
        raise ValueError("target_cardinality[j] - 1 < 0, ILP inválido.")

    # sparse_vstack devuelve un union type según el stub; pasamos por Any
    # antes de invocar tocsr() para evitar que Pylance pierda el método.
    a_stacked_any: Any = sparse_vstack([a_point, a_card])
    a_full = cast(csr_matrix, a_stacked_any.tocsr())
    b_full = np.concatenate([b_point, b_card])
    return a_full, b_full


# ============================================================================
# Algoritmo CSCLP completo
# ============================================================================

def run_csclp_algo(
    dataset: pd.DataFrame,
    target_cardinality: list[int],
    dist_method: str = "cosine",
    seed: int = 123,
) -> dict[str, Any]:
    """
    Ejecuta CSCLP sobre un dataset:
      1) prepara X/y
      2) calcula distancias N×N (se reutilizan después para el silhouette)
      3) PAM para encontrar k medoids
      4) ILP para asignar los no-medoids respetando cardinalidad
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

    # Validaciones (equivalentes a las del R original)
    if k < 2:
        raise ValueError("k inválido (<2).")
    if k >= n:
        raise ValueError(f"k={k} inválido para n={n}.")
    if np.any(np.isnan(target.astype(float))):
        raise ValueError("target_cardinality contiene NaN.")
    if int(target.sum()) != n:
        raise ValueError(
            f"target_cardinality suma {int(target.sum())} pero n={n}."
        )
    if np.any(target < 1):
        raise ValueError("target_cardinality tiene valores < 1.")

    # ----- Distancia completa N×N (también se reutiliza para silhouette) -----
    # Bug del stub de scipy: pdist sí acepta strings como `metric`. Casteamos
    # el retorno a ndarray para que squareform reciba un tipo conocido.
    condensed = cast(np.ndarray, pdist(x_mat, metric=dist_method))  # pyright: ignore[reportArgumentType, reportCallIssue]
    full_dist = squareform(condensed)

    # ----- PAM para k medoids -----
    # Llamamos a default_rng(seed) por consistencia con el set.seed(seed) del R
    # original (cluster::pam con random.start=FALSE no usa el RNG, pero
    # mantenemos el patrón para reproducibilidad en futuros cambios).
    _ = np.random.default_rng(seed)
    medoid_idx = pam_kaufman_rousseeuw(full_dist, k)

    # No-medoids
    all_idx        = np.arange(n)
    non_medoid_idx = np.setdiff1d(all_idx, medoid_idx, assume_unique=True)
    m              = int(len(non_medoid_idx))

    # ----- Vector de costos para el ILP -----
    # cost_mat[j, i] = distancia del medoid j al no-medoid i  (k × m)
    cost_mat = full_dist[np.ix_(medoid_idx, non_medoid_idx)]
    # Aplanado row-major: medoid 0 a todos, medoid 1 a todos, ...
    f_obj = cost_mat.flatten("C")

    # ----- Restricciones (sparse) -----
    a_full, b_full = _build_ilp_constraints(k, m, target)

    # ----- Resolver ILP con scipy.optimize.milp -----
    # Bug del stub: LinearConstraint y Bounds aceptan arrays para lb/ub,
    # pero el stub declara solo `float`. Silenciamos en la misma línea.
    constraints = LinearConstraint(a_full, lb=b_full, ub=b_full)  # pyright: ignore[reportArgumentType]
    integrality = np.ones(k * m, dtype=int)  # 1 = entero; con bounds [0,1] -> binario
    bounds      = Bounds(lb=0.0, ub=1.0)

    start_total = time.perf_counter()
    res = milp(
        c           = f_obj,
        constraints = constraints,
        integrality = integrality,
        bounds      = bounds,
    )
    total_time = time.perf_counter() - start_total

    # cast a Any porque OptimizeResult no tipa bien sus atributos dinámicos
    res_any = cast(Any, res)
    if not res_any.success or res_any.x is None:
        raise RuntimeError(
            f"ILP infeasible/failed. status={res_any.status} "
            f"message={res_any.message}"
        )

    # ----- Extraer solución -----
    x_sol = np.asarray(res_any.x, dtype=float)
    # Reshape a k × m: cada fila j es la asignación binaria del cluster j
    assign_mat = x_sol.reshape(k, m)

    y_pred = np.zeros(n, dtype=int)
    y_pred[medoid_idx] = np.arange(1, k + 1)        # medoids reciben labels 1..k
    # Equivalente a max.col(t(assign_mat)) de R: argmax por columna
    y_pred[non_medoid_idx] = np.argmax(assign_mat, axis=0) + 1

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
            "lp_solver":    "scipy.optimize.milp (HiGHS)",
            "all_bin":      True,
        },
    }


# ============================================================================
# Runner por dataset
# ============================================================================

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


# ============================================================================
# Punto de entrada para testing.py
# ============================================================================

def run(odatasets_unique: pd.DataFrame) -> None:
    """Loop principal equivalente al de CSCLP.R."""
    results: list[dict[str, Any]] = []

    for i in range(len(odatasets_unique)):
        print(
            f"\n--- Running CSCLP for dataset at position: {i + 1} ---",
            flush=True,
        )
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
                print(f"Time (LP): {row['Execution_Time']:.4f} s")
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
