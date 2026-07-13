"""
K-MeansACO
=================
ACO (Ant Colony Optimization).

Estructura del algoritmo:
  - 50 hormigas, 20 iteraciones
  - Matriz de feromonas tau de tamaño N x K (punto i -> cluster j)
  - La solución de K-Means se siembra como mejor candidata inicial: garantiza
    que ACO nunca devuelva un resultado peor que K-Means (ACO híbrido/sembrado)
  - Heurística eta[i, j] = 1 / (dist_coseno(punto i, centroide j) + eps)
  - Cada hormiga construye una asignación punto por punto eligiendo cluster con
    probabilidad proporcional a  (tau ** alpha) * (eta ** beta), respetando la
    capacidad de cada cluster (cardinalidad fija) -> soluciones siempre factibles
  - Función objetivo: silhouette (coseno) - penalty * |violación de cardinalidad|
    (con construcción capacitada la violación es 0, así que el objetivo es la
     silueta; la penalización se conserva como salvaguarda)
  - Actualización de feromonas: evaporación + depósito tipo "rank-based AS"
    (las mejores hormigas de la iteración y la mejor global depositan feromona)
  - Al final aplica adjust_cardinality a la mejor solución encontrada

"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
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
    get_predictions_dir,
    prepare_data,
    tabulate_clusters,
    to_int_list,
)

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(message)s")


# ============================================================================
# Hiperparámetros
# ============================================================================

@dataclass
class AcoHyperparams:
    """Hiperparámetros globales del algoritmo ACO."""
    algorithm:                   str   = "ACO"
    distance_method:             str   = "cosine"
    n_ants:                      int   = 50
    max_iterations:              int   = 20
    penalty_weight:              float = 100.0
    # --- feromonas ---
    alpha:                       float = 1.0     # peso de la feromona
    beta:                        float = 2.0     # peso de la heurística
    rho:                         float = 0.1     # tasa de evaporación
    n_rank:                      int   = 6       # nº de hormigas que depositan (rank-based)
    deposit_q:                   float = 1.0     # constante de depósito Q
    tau0:                        float = 1.0     # feromona inicial
    # --- siembra de la heurística (K-Means inicial) ---
    kmeans_nstart:               int   = 5
    kmeans_iter_max:             int   = 30
    adjust_cardinality_max_iter: int   = 1000
    master_seed:                 int   = 123
    objective:                   str   = "silhouette_cosine_minus_penalty"


ACO_HYPERPARAMS: dict[str, Any] = {
    "global": AcoHyperparams(),
    "runs":   {},
}

# ============================================================================
# Heurística: distancia coseno punto -> centroide
# ============================================================================

def _cosine_dist_to_centroids(
    x_mat: np.ndarray,
    centroids: np.ndarray,
    eps: float = 1e-12,
) -> np.ndarray:
    """
    Distancia coseno (N x K) de cada punto a cada centroide.

    """
    xn = x_mat / (np.linalg.norm(x_mat, axis=1, keepdims=True) + eps)
    cn = centroids / (np.linalg.norm(centroids, axis=1, keepdims=True) + eps)
    sim = xn @ cn.T                      # similitud coseno N x K
    return 1.0 - sim                     # distancia coseno N x K


def _heuristic_from_centroids(
    x_mat: np.ndarray,
    centroids: np.ndarray,
    eps: float = 1e-6,
) -> np.ndarray:
    """
    Heurística eta[i, j] = 1 / (dist_coseno(i, j) + eps).

    Cuanto más cerca está el punto i del centroide j, mayor es eta.
    Se normaliza por el máximo global para mantener magnitudes acotadas

    """
    dist_pc = _cosine_dist_to_centroids(x_mat, centroids)
    eta = 1.0 / (dist_pc + eps)
    max_eta = float(eta.max())
    if max_eta > 0.0 and np.isfinite(max_eta):
        eta = eta / max_eta
    return eta


# ============================================================================
# Construcción de solución guiada por feromonas (con capacidad por cluster)
# ============================================================================

def _construct_solution(
    desirability: np.ndarray,        # (tau ** alpha) * (eta ** beta), N x K
    target_cardinality: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Una hormiga construye una asignación punto por punto.

    Para cada punto (en orden aleatorio) elige un cluster con probabilidad
    proporcional a `desirability[i, :]`, pero solo entre clusters que aún
    tienen capacidad disponible. Como la suma de las
    cardinalidades objetivo es N (es la distribución real de clases), cada
    hormiga produce una solución EXACTAMENTE factible.

    Devuelve la asignación 1-indexada (valores en 1..k).
    """
    n, k = desirability.shape
    remaining = np.asarray(target_cardinality, dtype=int).copy()  # capacidad por cluster
    assignment = np.empty(n, dtype=int)

    order = rng.permutation(n)
    for i in order:
        probs = desirability[i].copy()
        full = remaining <= 0
        probs[full] = 0.0

        total = probs.sum()
        if total <= 0.0 or not np.isfinite(total):
            # Salvaguarda: si todo quedó en 0 (raro), asigna a un cluster con cupo
            avail = np.where(~full)[0]
            j = int(rng.choice(avail)) if avail.size else int(rng.integers(0, k))
        else:
            j = int(rng.choice(k, p=probs / total))

        assignment[i] = j + 1          # 1-indexado
        remaining[j] -= 1

    return assignment


# ============================================================================
# Actualización de feromonas: evaporación + depósito rank-based
# ============================================================================

def _deposit(tau: np.ndarray, ant: np.ndarray, amount: float) -> None:
    """Suma `amount` a las celdas tau[i, cluster(i)] usadas por la hormiga (in-place)."""
    rows = np.arange(len(ant))
    cols = ant - 1                       # de 1-indexado a 0-indexado
    np.add.at(tau, (rows, cols), amount)


def _score_to_quality(score: float) -> float:
    """
    Mapea el score (silueta en [-1, 1] cuando la solución es factible) a un
    factor de calidad positivo en [0, 1] para escalar el depósito.
    """
    if not np.isfinite(score):
        return 0.0
    return float(np.clip((score + 1.0) / 2.0, 0.0, 1.0))


def _update_pheromones(
    tau: np.ndarray,
    ranked: list[tuple[float, np.ndarray]],   # (score, ant) ordenadas DESC por score
    best_global: Optional[np.ndarray],
    best_global_score: float,
    rho: float,
    n_rank: int,
    deposit_q: float,
) -> None:
    """
    Actualización rank-based (AS_rank), in-place sobre tau:

      1) Evaporación:  tau <- (1 - rho) * tau
      2) Depósito:     las (n_rank - 1) mejores hormigas de la iteración
                       depositan con peso decreciente (n_rank-1, n_rank-2, ..., 1)
                       y la mejor solución global deposita con peso n_rank.
                       El depósito se escala por la calidad de cada solución.
    """
    tau *= (1.0 - rho)

    top = ranked[: max(0, n_rank - 1)]
    for r, (score, ant) in enumerate(top):
        weight = (n_rank - 1) - r        # r=0 -> n_rank-1 (la mejor de la iteración)
        amount = weight * deposit_q * _score_to_quality(score)
        if amount > 0.0:
            _deposit(tau, ant, amount)

    if best_global is not None:
        amount = n_rank * deposit_q * _score_to_quality(best_global_score)
        if amount > 0.0:
            _deposit(tau, best_global, amount)


# ============================================================================
# Centroides iniciales (siembra de la heurística vía K-Means)
# ============================================================================

def _initial_kmeans(
    x_mat: np.ndarray,
    k: int,
    rng: np.random.Generator,
    kmeans_nstart: int,
    kmeans_iter_max: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    K-Means una única vez. Devuelve:
      - centroids: K x D, para sembrar la heurística de la iteración 1
      - labels:    N, etiquetas 1-indexadas, para sembrar la mejor solución inicial

    """
    km = KMeans(
        n_clusters=k,
        n_init=kmeans_nstart,
        max_iter=kmeans_iter_max,
        random_state=int(rng.integers(0, 2**31 - 1)),
    )
    labels = km.fit_predict(x_mat).astype(int) + 1
    centroids = np.asarray(km.cluster_centers_, dtype=float)
    return centroids, labels


# ============================================================================
# ACO (loop principal con feromonas)
# ============================================================================

def run_aco_algorithm(
    x_mat: np.ndarray,
    target_cardinality: np.ndarray,
    n_ants:            int   = 50,
    max_iterations:    int   = 20,
    penalty_weight:    float = 100.0,
    alpha:             float = 1.0,
    beta:              float = 2.0,
    rho:               float = 0.1,
    n_rank:            int   = 6,
    deposit_q:         float = 1.0,
    tau0:              float = 1.0,
    kmeans_nstart:     int   = 5,
    kmeans_iter_max:   int   = 30,
    master_seed:       int   = 123,
) -> dict[str, Any]:
    """ACO con feromonas, heurística por centroides y construcción capacitada."""
    rng   = np.random.default_rng(master_seed)
    x_mat = np.asarray(x_mat, dtype=float)
    n     = int(x_mat.shape[0])
    k     = int(len(target_cardinality))
    target = np.asarray(target_cardinality, dtype=int)

    if k < 2:
        raise ValueError("k inválido (<2).")
    if k >= n:
        raise ValueError(f"k={k} inválido para n={n}.")

    # OPT: distancia coseno N×N una sola vez (para evaluar la silueta)
    d_cosine_sq = squareform(pdist(x_mat, metric="cosine"))

    # Feromonas inicializadas a un valor constante tau0
    tau = np.full((n, k), float(tau0), dtype=float)

    # Centroides + etiquetas iniciales de K-Means (una sola vez)
    centroids, km_labels = _initial_kmeans(
        x_mat, k, rng, kmeans_nstart, kmeans_iter_max
    )

    best_score:    float                = float("-inf")
    best_solution: Optional[np.ndarray] = None

    # --- SIEMBRA: la solución de K-Means como primera candidata ---------------
    # K-Means no respeta la cardinalidad, así que primero la hacemos factible
    # con adjust_cardinality y luego la evaluamos. Al fijarla como best_solution
    # inicial, ACO arranca desde una solución fuerte y queda garantizado a
    # terminar siendo >= K-Means (nunca peor). Además, esta solución se deposita
    # como "mejor global" en la primera actualización de feromonas, sesgando la
    # búsqueda hacia la región que K-Means encontró.
    seed_solution = adjust_cardinality(
        km_labels, x_mat, centroids, target_cardinality
    )
    seed_score = evaluate_solution(
        seed_solution, d_cosine_sq, target_cardinality, penalty_weight
    )
    if np.isfinite(seed_score):
        best_score    = float(seed_score)
        best_solution = np.asarray(seed_solution, dtype=int).copy()

    for _ in range(max_iterations):
        # 1) Heurística de esta iteración (a partir de los centroides actuales)
        eta = _heuristic_from_centroids(x_mat, centroids)

        # 2) Deseabilidad combinada (igual para todas las hormigas de la iter.)
        #    desirability[i, j] = tau[i, j]^alpha * eta[i, j]^beta
        desirability = np.power(tau, alpha) * np.power(eta, beta)

        # 3) Cada hormiga construye una solución factible guiada por feromonas
        scored: list[tuple[float, np.ndarray]] = []
        for _ in range(n_ants):
            ant = _construct_solution(desirability, target, rng)
            score = evaluate_solution(
                ant, d_cosine_sq, target_cardinality, penalty_weight
            )
            scored.append((score, ant))

            if np.isfinite(score) and score > best_score:
                best_score    = score
                best_solution = ant.copy()

        # 4) Ranking de hormigas de la iteración (mejor -> peor)
        scored.sort(key=lambda t: t[0], reverse=True)

        # 5) Actualización de feromonas (evaporación + depósito rank-based)
        _update_pheromones(
            tau, scored, best_solution, best_score,
            rho=rho, n_rank=n_rank, deposit_q=deposit_q,
        )

        # 6) Recalcular centroides desde la mejor solución global para que la
        #    heurística de la próxima iteración sea coherente con el óptimo actual
        if best_solution is not None:
            centroids = calculate_centroids(x_mat, best_solution, k)

    # Ajuste final de cardinalidad sobre la mejor solución (salvaguarda)
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
        "alpha":              1.0,
        "beta":               2.0,
        "rho":                0.1,
        "n_rank":             6,
        "deposit_q":          1.0,
        "tau0":               1.0,
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

        except Exception as e:
            print(f"Error: {e}")

    if results:
        df_out = pd.DataFrame(results)
        out_path = get_predictions_dir() / "pred_ACO.csv"
        df_out.to_csv(out_path, index=False)
        print(f"\nACO guardado en: {out_path}")
    else:
        print("\nACO finalizado sin resultados válidos. No se generó CSV.")
