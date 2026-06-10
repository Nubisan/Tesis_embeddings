"""
algorithms/sck1_final.py
========================
SCK1 (Size-Constrained K-means con ILP exacto).

Estructura del algoritmo:
  1. Inicialización: k centroides random (sample sin reemplazo de los puntos).
  2. Una pasada de ILP con cardinalidad EXACTA (igualdad):
     - Variables binarias x[i,j] = 1 si el punto j va al cluster i
     - Restr 1: cada punto a exactamente 1 cluster
     - Restr 2: cada cluster con exactamente sizes[i] puntos
     - Objetivo: minimizar suma de distancias coseno punto↔centroide
  3. Decodifica la asignación (1-based).
  4. (Calcula nuevos centroides pero no los usa, igual que el R original.)

Diferencias con CSCLP:
  - CSCLP optimiza solo sobre los m no-medoids (los medoids quedan fijos);
    SCK1 optimiza sobre todos los n puntos.
  - SCK1 inicializa con random sample, no con PAM.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, cast

import numpy as np
import pandas as pd
# scipy no provee stubs perfectos; silenciamos los warnings de tipos.
from scipy.optimize import Bounds, LinearConstraint, milp  # type: ignore[reportMissingTypeStubs]
from scipy.sparse import csr_matrix  # type: ignore[reportMissingTypeStubs]
from scipy.sparse import vstack as sparse_vstack  # type: ignore[reportMissingTypeStubs]
from scipy.spatial.distance import cdist, pdist, squareform

from ._common import (
    compute_metrics,
    get_predictions_dir,
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

@dataclass
class Sck1Hyperparams:
    """Hiperparámetros globales del algoritmo SCK1."""
    distance_method:   str = "cosine"
    ilp_direction:     str = "min"
    default_seed:      int = 123
    max_iterations:    int = 1
    centroid_init:     str = "random_sample"
    assignment_method: str = "ilp_exact_equality"
    centroid_update:   str = "colMeans"


SCK1_HYPERPARAMS: dict[str, Any] = {
    "global": Sck1Hyperparams(),
    "runs":   {},
}

# La ruta de salida se obtiene en runtime via get_predictions_dir()
# de _common, que lee la env var CLUSTERING_MODEL seteada por testing.py.


# ============================================================================
# Construcción del ILP
# ============================================================================

def _build_ilp_constraints(
    n: int, k: int, sizes: np.ndarray
) -> tuple[csr_matrix, np.ndarray]:
    """
    Construye las restricciones del ILP en formato sparse.

    Variables: x[i*n + j] = 1 si el punto j va al cluster i.
        i en [0, k), j en [0, n). Total: n_vars = k * n

    Indexación con orientación K × N (cluster × punto), igual que el R original
    (que transpone cost_mat antes de aplanarlo).

    Restricción 1 (n filas): cada punto a exactamente 1 cluster (igualdad)
        sum_i x[i*n + j] = 1   para todo j

    Restricción 2 (k filas): cada cluster con exactamente sizes[i] puntos
        sum_j x[i*n + j] = sizes[i]   para todo i

    Devuelve (A_full, b_full) donde lb=ub=b_full (todas igualdades).
    """
    n_vars = k * n

    # --- Restricción 1: A_point[j, i*n + j] = 1 ---
    # filas: rows_p[c] = c % n (índice de punto)
    # columnas: 0,1,...,n_vars-1
    rows_p = np.tile(np.arange(n), k)
    cols_p = np.arange(n_vars)
    data_p = np.ones(n_vars, dtype=float)
    a_point = csr_matrix((data_p, (rows_p, cols_p)), shape=(n, n_vars))
    b_point = np.ones(n, dtype=float)

    # --- Restricción 2: A_card[i, i*n + j] = 1 ---
    # filas: rows_c[c] = c // n (índice de cluster)
    rows_c = np.repeat(np.arange(k), n)
    cols_c = np.arange(n_vars)
    data_c = np.ones(n_vars, dtype=float)
    a_card = csr_matrix((data_c, (rows_c, cols_c)), shape=(k, n_vars))
    b_card = sizes.astype(float)

    # --- Combinar ambas restricciones ---
    a_stacked_any: Any = sparse_vstack([a_point, a_card])
    a_full = cast(csr_matrix, a_stacked_any.tocsr())
    b_full = np.concatenate([b_point, b_card])
    return a_full, b_full


def _solve_ilp_assignment(
    cost_mat_kn: np.ndarray,
    sizes: np.ndarray,
) -> np.ndarray:
    """
    Resuelve el ILP de asignación dada la matriz de costos K × N.

    cost_mat_kn[i, j] = costo de asignar el punto j al cluster i.

    Devuelve un vector p (1-based) con la asignación cluster de cada punto.
    """
    k = int(cost_mat_kn.shape[0])
    n = int(cost_mat_kn.shape[1])

    if len(sizes) != k:
        raise ValueError("sizes length != k")
    if int(sizes.sum()) != n:
        raise ValueError(
            f"sum(sizes)={int(sizes.sum())} != n={n} (inconsistente)"
        )
    if k < 2:
        raise ValueError("k < 2")

    # cost_vec aplanado row-major: cluster 0 a todos, cluster 1 a todos, ...
    # equivale a `as.vector(t(cost_mat))` de R con cost_mat K×N
    cost_vec = cost_mat_kn.flatten("C")

    # Construir restricciones (igualdades)
    a_full, b_full = _build_ilp_constraints(n, k, sizes)

    # Bug del stub: LinearConstraint acepta arrays para lb/ub.
    constraints = LinearConstraint(a_full, lb=b_full, ub=b_full)  # pyright: ignore[reportArgumentType]
    integrality = np.ones(k * n, dtype=int)
    bounds      = Bounds(lb=0.0, ub=1.0)

    res = milp(
        c           = cost_vec,
        constraints = constraints,
        integrality = integrality,
        bounds      = bounds,
    )

    # cast a Any porque OptimizeResult no tipa bien sus atributos dinámicos
    res_any = cast(Any, res)
    if not res_any.success or res_any.x is None:
        raise RuntimeError(
            f"ILP failed: status={res_any.status}, message={res_any.message}"
        )

    # Reshape a (k, n): fila i son las indicadoras del cluster i
    # equivale a `matrix(res$solution, nrow=k, byrow=TRUE)` de R
    sol = np.asarray(res_any.x, dtype=float).reshape(k, n)

    # Para cada punto j, asignación = argmax sobre las filas (clusters)
    # equivale a `max.col(t(assign_mat), "first")` de R
    p = np.argmax(sol, axis=0) + 1  # 1-based como R
    return p.astype(int)


# ============================================================================
# SCK1 iterativo (en realidad 1 sola pasada, igual que R)
# ============================================================================

def sck1_iterativo(
    data: np.ndarray,
    k: int,
    initial_centroids: np.ndarray,
    sizes: np.ndarray,
    dist_method: str = "cosine",
) -> dict[str, Any]:
    """
    Una pasada de SCK1: ILP + actualización de centroides (no usada).
    Devuelve dict con `assignments` (1-based) y `centroids` actualizados.
    """
    data      = np.asarray(data, dtype=float)
    centroids = np.asarray(initial_centroids, dtype=float).copy()
    n         = int(data.shape[0])

    if n != int(sizes.sum()):
        raise ValueError("nrow(data) != sum(sizes)")
    if int(centroids.shape[0]) != k:
        raise ValueError("nrow(centroids) != k")

    # Matriz de costos K × N (transpuesta de cdist N × K)
    # Bug del stub de scipy: cdist sí acepta strings como `metric`.
    cost_mat_nk = cast(np.ndarray, cdist(data, centroids, metric=dist_method))  # pyright: ignore[reportArgumentType, reportCallIssue]
    cost_mat_kn = cost_mat_nk.T

    assignments = _solve_ilp_assignment(cost_mat_kn, sizes)

    # Actualización de centroides (replicada por fidelidad con R, aunque
    # con max_iterations=1 los nuevos centroides nunca se usan).
    for i in range(1, k + 1):
        idx = np.where(assignments == i)[0]
        if len(idx) > 0:
            centroids[i - 1] = data[idx].mean(axis=0)

    return {"assignments": assignments, "centroids": centroids}


# ============================================================================
# Runner por dataset
# ============================================================================

def run_clustering_row(
    dataset: pd.DataFrame,
    target_cardinality: list[int],
    dataset_name: str,
    seed: int = 123,
) -> Optional[dict[str, Any]]:
    """Ejecuta SCK1 sobre un dataset y devuelve la fila de resultados."""
    data = prepare_data(dataset)
    if data is None:
        return None

    x_df:     pd.DataFrame    = data["X"]
    y_factor: pd.Categorical = data["y"]

    target = np.asarray(target_cardinality, dtype=int)
    x_mat  = x_df.to_numpy(dtype=float)
    n      = int(x_mat.shape[0])
    k      = int(len(target))

    if k < 2:
        return None
    if int(target.sum()) != n:
        return None

    start_total = time.perf_counter()

    # Inicialización: k centroides random (sample sin reemplazo de los puntos)
    rng = np.random.default_rng(seed)
    init_indices      = rng.choice(n, size=k, replace=False)
    initial_centroids = x_mat[init_indices].copy()

    try:
        result = sck1_iterativo(x_mat, k, initial_centroids, target)
    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.warning("Error en SCK1 para %s: %s", dataset_name, e)
        return None

    if result.get("assignments") is None:
        return None

    y_pred:   np.ndarray = result["assignments"]
    y_int                = (y_factor.codes + 1).astype(int)
    total_time           = time.perf_counter() - start_total

    # Pre-cálculo de la matriz de distancia N×N (reusada para silhouette,
    # igual que el R original con su comentario "OPT: distancia NxN una sola vez")
    # Bug del stub de scipy: pdist sí acepta strings como `metric`.
    condensed = cast(np.ndarray, pdist(x_mat, metric="cosine"))  # pyright: ignore[reportArgumentType, reportCallIssue]
    full_dist = squareform(condensed)

    metrics = compute_metrics(y_int, y_pred, x_mat, dist_matrix=full_dist)

    # Guardar hiperparámetros de la corrida
    SCK1_HYPERPARAMS["runs"][dataset_name] = {
        "seed":              int(seed),
        "k":                 k,
        "n":                 n,
        "n_features":        x_mat.shape[1],
        "target_sizes":      target.tolist(),
        "init_indices":      init_indices.tolist(),
        "distance_method":   "cosine",
        "max_iterations":    1,
        "centroid_init":     "random_sample",
        "assignment_method": "ilp_exact_equality",
    }

    cardinality_pred = tabulate_clusters(y_pred, k)

    return {
        "name":               dataset_name,
        "n":                  len(y_pred),
        "k":                  int(len(np.unique(y_pred))),
        "y_predict":          " ".join(str(int(c)) for c in y_pred),
        "y_true":             " ".join(str(int(c)) for c in y_int),
        "target_cardinality": " ".join(str(int(c)) for c in target),
        "cardinality_pred":   " ".join(str(int(c)) for c in cardinality_pred),
        "Execution_Time":     float(total_time),
        "ARI":                metrics["ARI"],
        "AMI":                metrics["AMI"],
        "NMI":                metrics["NMI"],
        "Silhouette_mean":    metrics["Silhouette_mean"],
    }


# ============================================================================
# Punto de entrada para testing.py
# ============================================================================

def run(odatasets_unique: pd.DataFrame) -> None:
    """Loop principal equivalente al de SCK1_final.R."""
    results: list[dict[str, Any]] = []

    for i in range(len(odatasets_unique)):
        print(
            f"\n--- Executing SCK1 for dataset at position: {i + 1} ---",
            flush=True,
        )
        try:
            dataset      = odatasets_unique["dataset"].iloc[i]
            dataset_name = odatasets_unique["name"].iloc[i]
            target_raw   = odatasets_unique["class_distribution_vector"].iloc[i]

            target = to_int_list(target_raw)
            if dataset is None or target is None or len(target) < 2:
                print("Dataset o cardinalidad inválida. Skipping.")
                continue

            row = run_clustering_row(dataset, target, dataset_name, seed=123)

            if row is not None:
                results.append(row)
                print(f"Time: {row['Execution_Time']:.4f} seconds")
            else:
                print("Dataset skipped.")

        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error: {e}")

    if results:
        df_out = pd.DataFrame(results)
        out_path = get_predictions_dir() / "pred_SCK1.csv"
        df_out.to_csv(out_path, index=False)
        print(f"\nSCK1 finalizado. Archivo guardado en: {out_path}")
    else:
        print("\nSCK1 finalizado sin resultados válidos. No se generó CSV.")
