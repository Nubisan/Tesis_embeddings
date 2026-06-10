"""
algorithms/pso.py
=================
PSO (Particle Swarm Optimization) con restricciones de cardinalidad.

Estructura del algoritmo:
  - Cada partícula es un vector de N reales en [0, 1].
  - Codificación order-based: ordenamos los valores y partimos en k segmentos
    de tamaños target[1], target[2], ..., target[k]. Garantiza cardinalidad
    exacta sin ajustes post-hoc (a diferencia de BAT y ACO).
  - Función de costo: suma de distancias intra-cluster (coseno).
  - 40 partículas × 20 iteraciones.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional, cast

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

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
class PsoHyperparams:
    """Hiperparámetros globales del algoritmo PSO."""
    algorithm:         str   = "PSO"
    distance_method:   str   = "cosine"
    pso_maxit:         int   = 20
    pso_swarm_size:    int   = 40
    par_lower:         float = 0.0
    par_upper:         float = 1.0
    objective:         str   = "intra_cluster_cosine_distance"
    assignment_method: str   = "order_based_partition"
    master_seed:       int   = 123


PSO_HYPERPARAMS: dict[str, Any] = {
    "global": PsoHyperparams(),
    "runs":   {},
}

# La ruta de salida se obtiene en runtime via get_predictions_dir()
# de _common, que lee la env var CLUSTERING_MODEL seteada por testing.py.


# ============================================================================
# Codificación order-based: vector real -> asignación de cluster
# ============================================================================

def _decode_par_to_clusters(
    par:             np.ndarray,
    k:               int,
    start_positions: np.ndarray,
    cum_sizes:       np.ndarray,
) -> np.ndarray:
    """
    Decodifica un vector real en [0,1] a una asignación de cluster (1-based).

    Ordena `par` ascendentemente y parte el vector ordenado en k segmentos
    de tamaños target[1], target[2], ..., target[k]:
      - posiciones [0, target[0])         -> cluster 1
      - posiciones [target[0], target[0]+target[1]) -> cluster 2
      - ...
    """
    cluster_assignment = np.zeros(len(par), dtype=int)
    ord_idx = np.argsort(par)
    for i in range(k):
        cluster_assignment[ord_idx[start_positions[i]:cum_sizes[i]]] = i + 1
    return cluster_assignment


# ============================================================================
# Función de costo: suma de distancias intra-cluster
# ============================================================================

def _cost_fn(
    par:             np.ndarray,
    dist_matrix:     np.ndarray,
    k:               int,
    start_positions: np.ndarray,
    cum_sizes:       np.ndarray,
) -> float:
    """
    Suma de distancias intra-cluster.

    Trick: como la matriz de distancias es simétrica con diagonal 0,
    sum(lower_tri) = sum(matriz) / 2. Esto es ~2x más rápido que np.tril.
    """
    cluster_assignment = _decode_par_to_clusters(par, k, start_positions, cum_sizes)
    total = 0.0
    for i in range(1, k + 1):
        ci = np.where(cluster_assignment == i)[0]
        if len(ci) > 1:
            cd = dist_matrix[np.ix_(ci, ci)]
            # cd simétrica con diagonal=0  =>  lower_tri.sum() = cd.sum() / 2
            total += float(cd.sum() * 0.5)
    return total


# ============================================================================
# PSO global-best (defaults SPSO 2007 de R psoptim)
# ============================================================================

def _pso_optimize(
    cost_fn:    Callable[[np.ndarray], float],
    n_dim:      int,
    swarm_size: int,
    max_iter:   int,
    lower:      float,
    upper:      float,
    rng:        np.random.Generator,
) -> tuple[np.ndarray, float]:
    """
    PSO global-best básico con defaults de R psoptim (SPSO 2007).

    Devuelve (mejor_posicion, mejor_score). Score se MINIMIZA (igual que R).
    """
    # Hiperparámetros estándar SPSO 2007 (defaults de R psoptim)
    w  = 1.0 / (2.0 * np.log(2.0))   # ≈ 0.7213
    c1 = 0.5 + np.log(2.0)            # ≈ 1.1931
    c2 = 0.5 + np.log(2.0)            # ≈ 1.1931

    # Inicializar posiciones uniformes en [lower, upper]
    positions:  np.ndarray = rng.uniform(lower, upper, size=(swarm_size, n_dim))
    velocities: np.ndarray = np.zeros((swarm_size, n_dim), dtype=float)

    # Evaluar posiciones iniciales
    pbest_scores: np.ndarray = np.array(
        [cost_fn(positions[i]) for i in range(swarm_size)],
        dtype=float,
    )
    pbest: np.ndarray = positions.copy()

    g_idx       = int(np.argmin(pbest_scores))
    gbest:      np.ndarray = pbest[g_idx].copy()
    gbest_score = float(pbest_scores[g_idx])

    for _ in range(max_iter):
        # Actualización de velocidades VECTORIZADA sobre todo el enjambre
        r1 = rng.random((swarm_size, n_dim))
        r2 = rng.random((swarm_size, n_dim))
        velocities = (
            w * velocities
            + c1 * r1 * (pbest - positions)
            + c2 * r2 * (gbest[np.newaxis, :] - positions)
        )
        positions = positions + velocities
        positions = np.clip(positions, lower, upper)

        # Evaluación: secuencial (cost_fn no se vectoriza fácilmente porque
        # cada partícula tiene un order distinto)
        for i in range(swarm_size):
            score = cost_fn(positions[i])
            if score < pbest_scores[i]:
                pbest[i]        = positions[i].copy()
                pbest_scores[i] = score
                if score < gbest_score:
                    gbest_score = score
                    gbest       = positions[i].copy()

    return gbest, gbest_score


# ============================================================================
# Algoritmo PSO completo
# ============================================================================

def run_pso_algorithm(
    x_mat: np.ndarray,
    target_cardinality: np.ndarray,
    swarm_size:  int   = 40,
    max_iter:    int   = 20,
    par_lower:   float = 0.0,
    par_upper:   float = 1.0,
    master_seed: int   = 123,
) -> dict[str, Any]:
    """
    Ejecuta PSO sobre x_mat. Devuelve dict con y_predict, best_score y full_dist.
    """
    rng    = np.random.default_rng(master_seed)
    x_mat  = np.asarray(x_mat, dtype=float)
    n      = int(x_mat.shape[0])
    target = np.asarray(target_cardinality, dtype=int)
    k      = int(len(target))

    # Validaciones equivalentes al R original
    if not np.all(np.isfinite(x_mat)):
        raise ValueError("X contiene NA/Inf")
    if k < 2:
        raise ValueError("k < 2")
    if int(target.sum()) != n:
        raise ValueError(
            f"sum(target_cardinality)={int(target.sum())} != n={n} (inconsistente)"
        )

    # ----- Distancia coseno N×N (reutilizada para silhouette luego) -----
    # Bug del stub de scipy: pdist sí acepta strings como `metric`.
    condensed = cast(np.ndarray, pdist(x_mat, metric="cosine"))  # pyright: ignore[reportArgumentType, reportCallIssue]
    full_dist = squareform(condensed)

    # ----- Pre-cálculo de cortes (OPT del R original, mantenido) -----
    cum_sizes       = np.cumsum(target)
    start_positions = np.concatenate([[0], cum_sizes[:-1]])

    # Closure que captura D, target, start_positions, cum_sizes
    def cost_fn(par: np.ndarray) -> float:
        return _cost_fn(par, full_dist, k, start_positions, cum_sizes)

    # ----- Optimización PSO -----
    best_par, best_score = _pso_optimize(
        cost_fn    = cost_fn,
        n_dim      = n,
        swarm_size = swarm_size,
        max_iter   = max_iter,
        lower      = par_lower,
        upper      = par_upper,
        rng        = rng,
    )

    # ----- Decodificar la mejor partícula -----
    y_pred = _decode_par_to_clusters(best_par, k, start_positions, cum_sizes)

    return {
        "y_predict":  y_pred,
        "best_score": best_score,
        "full_dist":  full_dist,
    }


# ============================================================================
# Runner por dataset
# ============================================================================

def run_clustering_row(
    dataset: pd.DataFrame,
    target_cardinality: list[int],
    dataset_name: str,
) -> Optional[dict[str, Any]]:
    """Ejecuta PSO sobre un dataset y devuelve la fila de resultados."""
    data = prepare_data(dataset)
    if data is None:
        return None

    x_df:     pd.DataFrame    = data["X"]
    y_factor: pd.Categorical = data["y"]

    target = np.asarray(target_cardinality, dtype=int)
    x_mat  = x_df.to_numpy(dtype=float)
    n      = int(x_mat.shape[0])
    k      = int(len(target))

    if k < 2 or int(target.sum()) != n:
        return None

    start_total = time.perf_counter()

    try:
        result = run_pso_algorithm(x_mat, target)
    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.warning("Error en run_pso_algorithm para %s: %s", dataset_name, e)
        return None

    y_pred:    np.ndarray = result["y_predict"]
    full_dist: np.ndarray = result["full_dist"]
    y_int                  = (y_factor.codes + 1).astype(int)
    total_time             = time.perf_counter() - start_total

    metrics = compute_metrics(y_int, y_pred, x_mat, dist_matrix=full_dist)

    PSO_HYPERPARAMS["runs"][dataset_name] = {
        "k":                 k,
        "n":                 n,
        "n_features":        x_mat.shape[1],
        "target_sizes":      target.tolist(),
        "pso_maxit":         20,
        "pso_swarm_size":    40,
        "par_lower":         0.0,
        "par_upper":         1.0,
        "distance_method":   "cosine",
        "assignment_method": "order_based_partition",
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
    """Loop principal equivalente al de PSO.R."""
    results: list[dict[str, Any]] = []

    for i in range(len(odatasets_unique)):
        print(
            f"\n--- Executing PSO for dataset at position: {i + 1} ---",
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
        out_path = get_predictions_dir() / "pred_PSO.csv"
        df_out.to_csv(out_path, index=False)
        print(f"\nPSO finalizado. Archivo guardado en: {out_path}")
    else:
        print("\nPSO finalizado sin resultados válidos. No se generó CSV.")
