"""
algorithms/km_milp.py
=====================
KM-MILP: K-Means con asignación vía Programación Lineal Entera (MILP)
con restricción de tamaño máximo por cluster.

Estructura del algoritmo:
  1. Inicialización: k centroides random (sample sin reemplazo de los puntos).
  2. Loop iterativo (máx 10 iter) hasta convergencia:
     - Asignación vía MILP: minimiza suma de distancias punto↔centroide,
       sujeto a (a) cada punto en exactamente 1 cluster, (b) tamaño por
       cluster <= size_constraints[j].
     - Actualización: nuevos centroides = media de cada cluster.
     - Convergencia: max(|c_old - c_new|) < 1e-4.

Diferencias con CSCLP:
  - CSCLP usa cardinalidad EXACTA (=); KM-MILP usa COTA SUPERIOR (<=).
  - CSCLP usa PAM con medoids fijos; KM-MILP itera centroides como kmeans.
"""

from __future__ import annotations

import logging
import os
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
class KmMilpHyperparams:
    """Hiperparámetros globales del algoritmo KM-MILP."""
    algorithm:        str   = "KM-MILP"
    distance_method:  str   = "cosine"
    seed:             int   = 123
    max_iter:         int   = 10
    convergence_tol:  float = 1e-4
    lp_solver:        str   = "scipy.optimize.milp (HiGHS)"
    all_bin:          bool  = True
    constraint_type:  str   = "<= (upper bound)"


KM_MILP_HYPERPARAMS: dict[str, Any] = {
    "global": KmMilpHyperparams(),
    "runs":   {},
}

# La ruta de salida se obtiene en runtime via get_predictions_dir()
# de _common, que lee la env var CLUSTERING_MODEL seteada por testing.py.


# ============================================================================
# Construcción del MILP
# ============================================================================

def _build_milp_constraints(
    n: int, k: int, size_constraints: np.ndarray
) -> tuple[csr_matrix, np.ndarray, np.ndarray]:
    """
    Construye las restricciones del MILP en formato sparse.

    Variables: x[i*k + j] = 1 si el punto i va al cluster j.
        i en [0, n), j en [0, k). Total: n_vars = n * k

    Restricción 1 (n filas): cada punto a exactamente 1 cluster (igualdad)
        sum_j x[i*k + j] = 1   para todo i

    Restricción 2 (k filas): cada cluster a lo sumo size[j] (cota superior)
        sum_i x[i*k + j] <= size[j]   para todo j

    Devuelve (A_full, lb_full, ub_full) donde
      lb=ub=1     para las primeras n filas (igualdad)
      lb=-inf, ub=size  para las últimas k filas (<=)
    """
    n_vars = n * k

    # --- Restricción 1: A_point[i, i*k + j] = 1 ---
    # filas: cada índice de punto se repite k veces
    # columnas: 0,1,...,n*k-1
    rows_p = np.repeat(np.arange(n), k)
    cols_p = np.arange(n_vars)
    data_p = np.ones(n_vars, dtype=float)
    a_point = csr_matrix((data_p, (rows_p, cols_p)), shape=(n, n_vars))

    # --- Restricción 2: A_card[j, i*k + j] = 1 (todas las filas con ese j) ---
    # filas: 0,1,...,k-1 se repiten en bloques de tamaño k (uno por punto)
    rows_c = np.tile(np.arange(k), n)
    cols_c = np.arange(n_vars)
    data_c = np.ones(n_vars, dtype=float)
    a_card = csr_matrix((data_c, (rows_c, cols_c)), shape=(k, n_vars))

    # --- Combinar ambas restricciones ---
    a_stacked_any: Any = sparse_vstack([a_point, a_card])
    a_full = cast(csr_matrix, a_stacked_any.tocsr())

    # bounds: igualdad para las primeras n filas, cota superior para las últimas k
    lb_full = np.concatenate([np.ones(n), np.full(k, -np.inf)])
    ub_full = np.concatenate([np.ones(n), size_constraints.astype(float)])

    return a_full, lb_full, ub_full


def _solve_milp_assignment(
    data: np.ndarray,
    centroids: np.ndarray,
    size_constraints: np.ndarray,
    a_full: csr_matrix,
    lb_full: np.ndarray,
    ub_full: np.ndarray,
    integrality: np.ndarray,
    bounds: Bounds,
    dist_method: str,
    time_limit: Optional[float] = None,
) -> tuple[np.ndarray, str]:
    """
    Resuelve el MILP de asignación dados los centroides actuales.
    Devuelve (p, status) con p = vector 1-based de asignación y status en
    {"optimal", "time_limit"}.

    Si `time_limit` (segundos) se pasa y el solver lo alcanza sin probar
    optimalidad, HiGHS devuelve la MEJOR solución factible encontrada
    (incumbente): se acepta y se marca status="time_limit". Solo se lanza
    error si NO se encontró ninguna solución factible.

    Las matrices de restricciones se reciben pre-construidas para no
    rehacerlas en cada iteración del loop kmeans-style.
    """
    n = int(data.shape[0])
    k = int(centroids.shape[0])

    if k < 2 or n <= 1:
        raise ValueError(f"n={n}, k={k} inválidos.")
    if len(size_constraints) != k:
        raise ValueError("size_constraints length != k")

    # Costo: distancia de cada punto a cada centroide (cdist usa BLAS)
    # Bug del stub de scipy: cdist sí acepta strings como `metric`.
    cost_mat = cast(np.ndarray, cdist(data, centroids, metric=dist_method))  # pyright: ignore[reportArgumentType, reportCallIssue]
    # Aplanado row-major: punto 0 a clusters 0..k-1, punto 1 a clusters 0..k-1, ...
    cost_vec = cost_mat.flatten("C")

    # Bug del stub: LinearConstraint acepta arrays para lb/ub.
    constraints = LinearConstraint(a_full, lb=lb_full, ub=ub_full)  # pyright: ignore[reportArgumentType]

    options = {"time_limit": float(time_limit)} if time_limit else None
    res = milp(
        c           = cost_vec,
        constraints = constraints,
        integrality = integrality,
        bounds      = bounds,
        options     = options,
    )

    # cast a Any porque OptimizeResult no tipa bien sus atributos dinámicos
    res_any = cast(Any, res)
    # status 0 = óptimo; status 1 = límite (tiempo/iter) alcanzado. En ambos,
    # si hay solución factible (x no es None) la usamos.
    if res_any.x is None:
        raise RuntimeError(
            f"MILP sin solución factible: status={res_any.status}, "
            f"message={res_any.message}"
        )
    status = "optimal" if res_any.success else "time_limit"

    # Reshape a (n, k): cada fila i son las k indicadoras del punto i
    sol = np.asarray(res_any.x, dtype=float).reshape(n, k)
    # max.col(sol, "first") en R == argmax por fila en Python (con ties.first)
    p = np.argmax(sol, axis=1) + 1  # 1-based como R
    return p.astype(int), status


# ============================================================================
# Loop iterativo tipo KMeans con asignación MILP
# ============================================================================

def clustering_with_size_constraints(
    data: np.ndarray,
    size_constraints: np.ndarray,
    max_iter:        int   = 10,
    seed:            int   = 123,
    dist_method:     str   = "cosine",
    convergence_tol: float = 1e-4,
    time_limit:      Optional[float] = None,
) -> dict[str, Any]:
    """
    K-Means con asignación MILP. Itera hasta convergencia o max_iter.

    `time_limit` (segundos, por resolución MILP) es opcional: si se alcanza sin
    probar optimalidad, se usa la mejor solución factible (incumbente) y se
    marca la corrida como aproximada.

    Devuelve dict con `p` (asignación 1-based) e `hyperparams`.
    """
    data = np.asarray(data, dtype=float)
    n    = int(data.shape[0])
    k    = int(len(size_constraints))

    if k < 2:
        raise ValueError("k < 2")
    if int(size_constraints.sum()) < n:
        raise ValueError(
            f"sum(size_constraints)={int(size_constraints.sum())} < n={n}"
        )

    # Inicialización: k centroides aleatorios (sample sin reemplazo)
    rng = np.random.default_rng(seed)
    init_indices = rng.choice(n, size=k, replace=False)
    centroids: np.ndarray = data[init_indices].copy()

    # Pre-construir las restricciones del MILP (no cambian entre iteraciones)
    a_full, lb_full, ub_full = _build_milp_constraints(n, k, size_constraints)
    integrality              = np.ones(n * k, dtype=int)
    bounds                   = Bounds(lb=0.0, ub=1.0)

    p: np.ndarray = np.ones(n, dtype=int)
    converged_iter = max_iter
    any_time_limited = False

    for iter_num in range(1, max_iter + 1):
        p_new, solve_status = _solve_milp_assignment(
            data, centroids, size_constraints,
            a_full, lb_full, ub_full,
            integrality, bounds,
            dist_method,
            time_limit=time_limit,
        )
        if solve_status == "time_limit":
            any_time_limited = True

        # Actualización de centroides como media por cluster
        new_centroids = centroids.copy()
        for j in range(1, k + 1):
            idx = np.where(p_new == j)[0]
            if len(idx) > 0:
                new_centroids[j - 1] = data[idx].mean(axis=0)

        # Verificar convergencia: max diff entre centroides viejos y nuevos
        if float(np.max(np.abs(centroids - new_centroids))) < convergence_tol:
            converged_iter = iter_num
            p = p_new
            break

        centroids = new_centroids
        p = p_new

    return {
        "p": p,
        "hyperparams": {
            "seed":              int(seed),
            "dist_method":       dist_method,
            "max_iter":          max_iter,
            "converged_at_iter": converged_iter,
            "convergence_tol":   convergence_tol,
            "init_indices":      init_indices.tolist(),
            "lp_solver":         "scipy.optimize.milp (HiGHS)",
            "all_bin":           True,
            "constraint_type":   "<= (upper bound)",
            "time_limit_s":      float(time_limit) if time_limit else None,
            "solution_type":     "incumbent" if any_time_limited else "optimal",
        },
    }


# ============================================================================
# Runner por dataset
# ============================================================================

def _registrar_incumbente(name: str, n: int, k: int,
                          elapsed: float, time_limit: Optional[float]) -> None:
    """
    Deja constancia (sidecar CSV, no toca el esquema de pred_KM-MILP.csv) de que
    un dataset se resolvió con la incumbente por límite de tiempo. Sirve para tu
    nota de fidelidad: distingue resultados exactos de aproximados.
    """
    try:
        out = get_predictions_dir() / "km_milp_incumbentes.csv"
        row = pd.DataFrame([{
            "timestamp":     pd.Timestamp.now().isoformat(),
            "model":         os.environ.get("CLUSTERING_MODEL", ""),
            "modalidad":     os.environ.get("CLUSTERING_MODALIDAD", ""),
            "dataset":       name,
            "n":             n,
            "k":             k,
            "elapsed_s":     round(elapsed, 1),
            "time_limit_s":  time_limit,
            "solution_type": "incumbent",
        }])
        row.to_csv(out, mode="a", header=not out.exists(), index=False)
    except Exception:  # pylint: disable=broad-exception-caught
        pass


def run_clustering_row(
    dataset: pd.DataFrame,
    target_cardinality: list[int],
    dataset_name: str,
) -> Optional[dict[str, Any]]:
    """Ejecuta KM-MILP sobre un dataset y devuelve la fila de resultados."""
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

    # time_limit por resolución MILP: se toma de la env var KM_MILP_TIME_LIMIT
    # (segundos). Si no está, comportamiento original (sin límite = exacto).
    _tl_raw = os.environ.get("KM_MILP_TIME_LIMIT", "").strip()
    time_limit = float(_tl_raw) if _tl_raw else None

    try:
        result = clustering_with_size_constraints(x_mat, target, time_limit=time_limit)
    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.warning("Error en KM-MILP para %s: %s", dataset_name, e)
        return None

    y_pred:    np.ndarray = result["p"]
    hp:        dict[str, Any] = result["hyperparams"]
    y_int                  = (y_factor.codes + 1).astype(int)
    total_time             = time.perf_counter() - start_total

    # Si se usó la incumbente (límite de tiempo), lo dejamos registrado en un
    # sidecar (NO cambia el esquema del CSV) para tu nota de fidelidad: así sabes
    # exactamente qué resultados de KM-MILP son aproximados y cuáles exactos.
    if hp.get("solution_type") == "incumbent":
        _registrar_incumbente(dataset_name, n, k, total_time, time_limit)
        logger.warning("KM-MILP %s: solución INCUMBENTE (límite %.0fs), aproximada.",
                       dataset_name, time_limit or 0.0)

    # Pre-cálculo de la matriz de distancia N×N (reusada para silhouette,
    # igual que el R original con su comentario "reuse instead of recomputing")
    # Bug del stub de scipy: pdist sí acepta strings como `metric`.
    condensed = cast(np.ndarray, pdist(x_mat, metric=hp["dist_method"]))  # pyright: ignore[reportArgumentType, reportCallIssue]
    full_dist = squareform(condensed)

    metrics = compute_metrics(y_int, y_pred, x_mat, dist_matrix=full_dist)

    KM_MILP_HYPERPARAMS["runs"][dataset_name] = {
        **hp,
        "k":            k,
        "n":            n,
        "n_features":   x_mat.shape[1],
        "target_sizes": target.tolist(),
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
    """Loop principal equivalente al de KM-MILP.R."""
    results: list[dict[str, Any]] = []

    for i in range(len(odatasets_unique)):
        print(
            f"\n--- Executing KM-MILP for dataset at position: {i + 1} ---",
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

            row = run_clustering_row(dataset, target, dataset_name)

            if row is not None:
                results.append(row)
                print(f"Time: {row['Execution_Time']:.4f} seconds")
            else:
                print("Dataset skipped.")

        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error: {e}")

    if results:
        df_out = pd.DataFrame(results)
        out_path = get_predictions_dir() / "pred_KM-MILP.csv"
        df_out.to_csv(out_path, index=False)
        print(f"\nKM-MILP finalizado. Archivo guardado en: {out_path}")
    else:
        print("\nKM-MILP finalizado sin resultados válidos. No se generó CSV.")