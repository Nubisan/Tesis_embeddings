"""
algorithms/aco.py
=================
ACO (Ant Colony Optimization) para clustering con cardinalidad fija.
Migración Python del archivo ACO.R.

Estructura del algoritmo:
  - 50 hormigas, 20 iteraciones
  - Cada hormiga genera una solución NUEVA en cada iteración (kmeans + perturbación)
  - Función objetivo: silhouette (coseno) - penalty * |violación de cardinalidad|
    (penalty_weight = 100, mucho más agresivo que BAT que usa 10)
  - Al final aplica adjust_cardinality a la mejor solución encontrada

Diferencia con R:
  - El comentario en ACO.R dice: "Loop original preservado: runif(1) por punto
    mantiene la secuencia exacta del RNG" (no vectorizaron la perturbación
    para mantener fidelidad RNG entre versiones de R).
  - En Python ya asumimos que el RNG es distinto a R, así que VECTORIZAMOS
    la perturbación. Esto da una mejora notable de velocidad: en lugar de N
    iteraciones del bucle, una sola operación con máscara booleana.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans

from ._common import (
    adjust_cardinality,
    calculate_centroids,
    compute_metrics,
    evaluate_solution,
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
class AcoHyperparams:
    """Hiperparámetros globales del algoritmo ACO."""
    algorithm:                   str   = "ACO"
    distance_method:             str   = "cosine"
    n_ants:                      int   = 50
    max_iterations:              int   = 20
    penalty_weight:              float = 100.0
    perturbation_rate:           float = 0.1
    kmeans_nstart:               int   = 5
    kmeans_iter_max:             int   = 30
    adjust_cardinality_max_iter: int   = 1000
    master_seed:                 int   = 123
    objective:                   str   = "silhouette_cosine_minus_penalty"


ACO_HYPERPARAMS: dict[str, Any] = {
    "global": AcoHyperparams(),
    "runs":   {},
}

PROJECT_ROOT    = Path(__file__).resolve().parent.parent
PREDICTIONS_DIR = PROJECT_ROOT / "predictions"
PREDICTIONS_DIR.mkdir(parents=True, exist_ok=True)
PRED_ACO_CSV    = PREDICTIONS_DIR / "pred_ACO.csv"


# ============================================================================
# Solución inicial: K-means + ajuste random simple
# ============================================================================

def _generate_initial_solution(
    x_mat: np.ndarray,
    target_cardinality: np.ndarray,
    rng: np.random.Generator,
    kmeans_nstart: int = 5,
    kmeans_iter_max: int = 30,
) -> np.ndarray:
    """
    Solución inicial: K-means + ajuste heurístico simple.

    A diferencia de BAT (que usa centroides para mover el punto al cluster
    más cercano), aquí ACO usa ajuste random: si un cluster está sobre-poblado,
    toma un punto al azar y lo manda a un cluster random distinto.
    """
    k = int(len(target_cardinality))

    # KMeans con random_state derivado del RNG para reproducibilidad encadenada
    km = KMeans(
        n_clusters=k,
        n_init=kmeans_nstart,
        max_iter=kmeans_iter_max,
        random_state=int(rng.integers(0, 2**31 - 1)),
    )
    cluster_assignment: np.ndarray = km.fit_predict(x_mat).astype(int) + 1

    target = np.asarray(target_cardinality, dtype=int)
    other_clusters_cache = {
        j: np.array([c for c in range(1, k + 1) if c != j], dtype=int)
        for j in range(1, k + 1)
    }

    for j in range(1, k + 1):
        while int(np.sum(cluster_assignment == j)) > int(target[j - 1]):
            idx = np.where(cluster_assignment == j)[0]
            if len(idx) == 0:
                break
            chosen_point = int(rng.choice(idx))
            new_cluster  = int(rng.choice(other_clusters_cache[j]))
            cluster_assignment[chosen_point] = new_cluster

    return cluster_assignment


# ============================================================================
# Perturbación vectorizada (mejora vs R)
# ============================================================================

def _perturb_solution(
    cluster_assignment: np.ndarray,
    k: int,
    perturbation_rate: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Cambia ~`perturbation_rate` de los puntos a un cluster aleatorio.

    R usa un loop punto por punto con runif(1) para preservar la secuencia
    del RNG. En Python vectorizamos: una máscara booleana + asignación masiva.
    """
    new_ca = cluster_assignment.copy()
    n      = int(len(new_ca))
    mask   = rng.random(n) < perturbation_rate
    if mask.any():
        new_ca[mask] = rng.integers(1, k + 1, size=int(mask.sum()))
    return new_ca


# ============================================================================
# ACO (loop principal)
# ============================================================================

def run_aco_algorithm(
    x_mat: np.ndarray,
    target_cardinality: np.ndarray,
    n_ants:            int   = 50,
    max_iterations:    int   = 20,
    penalty_weight:    float = 100.0,
    perturbation_rate: float = 0.1,
    kmeans_nstart:     int   = 5,
    kmeans_iter_max:   int   = 30,
    master_seed:       int   = 123,
) -> dict[str, Any]:
    """ACO con caché de distancia coseno y perturbación vectorizada."""
    rng   = np.random.default_rng(master_seed)
    x_mat = np.asarray(x_mat, dtype=float)
    n     = int(x_mat.shape[0])
    k     = int(len(target_cardinality))

    if k < 2:
        raise ValueError("k inválido (<2).")
    if k >= n:
        raise ValueError(f"k={k} inválido para n={n}.")

    # OPT: distancia coseno N×N una sola vez (ya estaba en R)
    d_cosine_sq = squareform(pdist(x_mat, metric="cosine"))

    best_score:    float                = float("-inf")
    best_solution: Optional[np.ndarray] = None

    for _ in range(max_iterations):
        # En cada iteración se genera una población NUEVA (sin memoria
        # entre iteraciones, a diferencia de BAT). Por eso ACO es más caro.
        ants: list[np.ndarray] = []
        for _ in range(n_ants):
            ca = _generate_initial_solution(
                x_mat, target_cardinality, rng,
                kmeans_nstart=kmeans_nstart,
                kmeans_iter_max=kmeans_iter_max,
            )
            ca = _perturb_solution(ca, k, perturbation_rate, rng)
            ants.append(ca)

        for ant in ants:
            score = evaluate_solution(
                ant, d_cosine_sq, target_cardinality, penalty_weight
            )
            if np.isfinite(score) and score > best_score:
                best_score    = score
                best_solution = ant.copy()

    # Ajuste final de cardinalidad sobre la mejor solución (igual que R)
    if best_solution is not None:
        centroids     = calculate_centroids(x_mat, best_solution, k)
        best_solution = adjust_cardinality(
            best_solution, x_mat, centroids, target_cardinality
        )

    return {
        "best_solution": best_solution,
        "best_score":    best_score,
    }


# ============================================================================
# Runner por dataset
# ============================================================================

def run_clustering_row(
    dataset: pd.DataFrame,
    target_cardinality: list[int],
    dataset_name: str,
) -> Optional[dict[str, Any]]:
    """Ejecuta ACO sobre un dataset y devuelve la fila de resultados."""
    data = prepare_data(dataset)
    if data is None:
        return None

    x_df:     pd.DataFrame    = data["X"]
    y_factor: pd.Categorical = data["y"]

    target = np.asarray(target_cardinality, dtype=int)
    x_mat  = x_df.to_numpy(dtype=float)

    start_total = time.perf_counter()

    try:
        result = run_aco_algorithm(x_mat, target)
    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.warning("Error en run_aco_algorithm para %s: %s", dataset_name, e)
        return None

    if result.get("best_solution") is None:
        return None

    y_pred     = np.asarray(result["best_solution"], dtype=int)
    y_int      = (y_factor.codes + 1).astype(int)
    total_time = time.perf_counter() - start_total

    metrics = compute_metrics(y_int, y_pred, x_mat)

    ACO_HYPERPARAMS["runs"][dataset_name] = {
        "k":                  len(target),
        "n":                  len(y_pred),
        "n_features":         x_mat.shape[1],
        "target_sizes":       target.tolist(),
        "n_ants":             50,
        "max_iterations":     20,
        "penalty_weight":     100.0,
        "perturbation_rate":  0.1,
        "distance_method":    "cosine",
    }

    cardinality_pred = tabulate_clusters(y_pred, len(target))

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
    """Loop principal equivalente al de ACO.R."""
    results: list[dict[str, Any]] = []

    for i in range(len(odatasets_unique)):
        print(
            f"\n--- Executing ACO for dataset at position: {i + 1} ---",
            flush=True,
        )
        try:
            dataset      = odatasets_unique["dataset"].iloc[i]
            dataset_name = odatasets_unique["name"].iloc[i]
            target_raw   = odatasets_unique["class_distribution_vector"].iloc[i]

            target = to_int_list(target_raw)
            if target is None or len(target) < 2:
                print("Target cardinality inválida. Skipping.")
                continue

            row = run_clustering_row(dataset, target, dataset_name)

            if row is not None:
                results.append(row)
                print(f"Time: {row['Execution_Time']:.4f} s")
            else:
                print("Dataset skipped.")

        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error: {e}")

    if results:
        df_out = pd.DataFrame(results)
        df_out.to_csv(PRED_ACO_CSV, index=False)
        print(f"\nACO guardado en: {PRED_ACO_CSV}")
    else:
        print("\nACO finalizado sin resultados válidos. No se generó CSV.")
