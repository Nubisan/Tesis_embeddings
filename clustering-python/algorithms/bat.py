"""
algorithms/bat.py
=================
Algoritmo de clustering BAT (metaheurístico inspirado en murciélagos).

Estructura del algoritmo:
  - 30 murciélagos exploran el espacio durante 20 iteraciones
  - Inicialización: K-means + ajuste de cardinalidad por centroide cercano
  - Función objetivo: silhouette (coseno) - penalty por desviación de cardinalidad
  - Distancia coseno precomputada N×N para evitar recálculos
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
    get_predictions_dir,
    prepare_data,
    tabulate_clusters,
    to_int_list,
)

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(message)s")


# ============================================================================
# Hiperparámetros y estado del murciélago
# ============================================================================

@dataclass
class BatHyperparams:
    """Hiperparámetros globales del algoritmo BAT."""
    algorithm:                   str   = "BAT"
    distance_method:             str   = "cosine"
    n_bats:                      int   = 30
    max_iterations:              int   = 20
    f_min:                       float = 0.0
    f_max:                       float = 2.0
    loudness_init:               float = 0.5
    pulse_rate_init:             float = 0.5
    alpha:                       float = 0.9
    gamma:                       float = 0.9
    penalty_weight:              float = 10.0
    master_seed:                 int   = 1521
    kmeans_nstart:               int   = 5
    kmeans_iter_max:             int   = 30
    adjust_cardinality_max_iter: int   = 1000
    objective:                   str   = "silhouette_cosine_minus_penalty"


@dataclass
class _BatState:
    """
    Estado mutable de un murciélago. Usamos dataclass (no TypedDict) porque
    Pylance propaga mejor los tipos a través de acceso por atributo que por
    clave de TypedDict.
    """
    position:   np.ndarray
    velocity:   np.ndarray
    frequency:  float
    loudness:   float
    pulse_rate: float
    seed:       int


BAT_HYPERPARAMS: dict[str, Any] = {
    "global": BatHyperparams(),
    "runs":   {},
}

# La ruta de salida se obtiene en runtime via get_predictions_dir()
# de _common, que lee la env var CLUSTERING_MODEL seteada por testing.py.


# ============================================================================
# Solución inicial: K-means + ajuste de cardinalidad
# ============================================================================

def generate_initial_solution(
    x_mat: np.ndarray,
    target_cardinality: np.ndarray,
    seed: int = 45,
    kmeans_nstart: int = 5,
    kmeans_iter_max: int = 30,
) -> np.ndarray:
    """K-means con `nstart` inicios + ajuste de cardinalidad por centroide cercano."""
    x_mat = np.asarray(x_mat, dtype=float)
    k     = len(target_cardinality)

    km = KMeans(
        n_clusters=k,
        n_init=kmeans_nstart,
        max_iter=kmeans_iter_max,
        random_state=int(seed),
    )
    labels = km.fit_predict(x_mat) + 1  # 1..k para alinearnos con R

    centroids = calculate_centroids(x_mat, labels, k)
    return adjust_cardinality(labels, x_mat, centroids, target_cardinality)


# ============================================================================
# BAT (loop principal)
# ============================================================================

def run_bat_algorithm(
    x_mat: np.ndarray,
    target_cardinality: np.ndarray,
    n_bats:         int   = 30,
    max_iterations: int   = 20,
    f_min:          float = 0.0,
    f_max:          float = 2.0,
    loudness:       float = 0.5,
    pulse_rate:     float = 0.5,
    alpha:          float = 0.9,
    gamma:          float = 0.9,
    penalty_weight: float = 10.0,
    master_seed:    int   = 1521,
) -> dict[str, Any]:
    """Algoritmo BAT con caché de distancia coseno y operaciones vectorizadas."""
    rng   = np.random.default_rng(master_seed)
    x_mat = np.asarray(x_mat, dtype=float)
    n     = int(x_mat.shape[0])
    k     = int(len(target_cardinality))

    # OPT: distancia coseno NxN se calcula una sola vez
    d_cosine_sq = squareform(pdist(x_mat, metric="cosine"))

    # Semillas únicas para los murciélagos (R: sample(1:10000, n_bats))
    seeds = rng.choice(np.arange(1, 10001), size=n_bats, replace=False)

    # Inicializar la población
    bats: list[_BatState] = []
    for i in range(n_bats):
        bats.append(_BatState(
            position   = generate_initial_solution(
                x_mat, target_cardinality, seed=int(seeds[i])
            ),
            velocity   = np.zeros(n, dtype=float),
            frequency  = float(rng.uniform(f_min, f_max)),
            loudness   = loudness,
            pulse_rate = pulse_rate,
            seed       = int(seeds[i]),
        ))

    best_solution: np.ndarray = bats[0].position.copy()
    best_score:    float       = evaluate_solution(
        best_solution, d_cosine_sq, target_cardinality, penalty_weight
    )
    best_seed:     int         = bats[0].seed

    for iteration in range(1, max_iterations + 1):
        for i in range(n_bats):
            bat = bats[i]
            bat.frequency = float(rng.uniform(f_min, f_max))
            bat.velocity  = np.asarray(
                bat.velocity + (bat.position - best_solution) * bat.frequency,
                dtype=float,
            )

            raw_pos                  = np.asarray(bat.position + bat.velocity, dtype=float)
            new_position: np.ndarray = np.asarray(np.round(raw_pos), dtype=int)
            new_position             = np.asarray(np.clip(new_position, 1, k), dtype=int)

            if rng.random() > bat.pulse_rate:
                new_position = np.asarray(rng.integers(1, k + 1, size=n), dtype=int)

            # Ajustar excesos por cardinalidad (lógica equivalente al while de R)
            for j in range(1, k + 1):
                while int(np.sum(new_position == j)) > int(target_cardinality[j - 1]):
                    idx = np.where(new_position == j)[0]
                    if len(idx) > 0:
                        new_position[int(rng.choice(idx))] = int(rng.integers(1, k + 1))
                    else:
                        break

            new_score = evaluate_solution(
                new_position, d_cosine_sq, target_cardinality, penalty_weight
            )

            if (
                np.isfinite(new_score)
                and new_score > best_score
                and rng.random() < bat.loudness
            ):
                bats[i].position   = new_position
                bats[i].loudness   = alpha * bat.loudness
                bats[i].pulse_rate = pulse_rate * (
                    1.0 - float(np.exp(-gamma * iteration))
                )

                best_solution = new_position.copy()
                best_score    = new_score
                best_seed     = bats[i].seed

    return {
        "best_solution": best_solution,
        "best_score":    best_score,
        "best_seed":     best_seed,
        "seeds":         seeds.tolist(),
    }


# ============================================================================
# Runner por dataset
# ============================================================================

def run_clustering_row(
    dataset: pd.DataFrame,
    target_cardinality: list[int],
    dataset_name: str,
) -> Optional[dict[str, Any]]:
    """Ejecuta BAT sobre un dataset y devuelve la fila de resultados."""
    data = prepare_data(dataset)
    if data is None:
        return None

    x_df:     pd.DataFrame    = data["X"]
    y_factor: pd.Categorical = data["y"]

    target = np.asarray(target_cardinality, dtype=int)
    x_mat  = x_df.to_numpy(dtype=float)

    start_total = time.perf_counter()

    try:
        result = run_bat_algorithm(x_mat, target)
    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.warning("Error en run_bat_algorithm para %s: %s", dataset_name, e)
        return None

    y_pred     = np.asarray(result["best_solution"], dtype=int)
    y_int      = (y_factor.codes + 1).astype(int)
    total_time = time.perf_counter() - start_total

    metrics = compute_metrics(y_int, y_pred, x_mat)

    BAT_HYPERPARAMS["runs"][dataset_name] = {
        "k":                len(target),
        "n":                len(y_pred),
        "n_features":       x_mat.shape[1],
        "target_sizes":     target.tolist(),
        "n_bats":           30,
        "max_iterations":   20,
        "f_min":            0.0,
        "f_max":            2.0,
        "loudness_init":    0.5,
        "pulse_rate_init":  0.5,
        "alpha":            0.9,
        "gamma":            0.9,
        "penalty_weight":   10.0,
        "master_seed":      1521,
        "best_seed":        result["best_seed"],
        "all_seeds":        result["seeds"],
        "distance_method":  "cosine",
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
    """Loop principal equivalente al de Bat.R."""
    results: list[dict[str, Any]] = []

    for i in range(len(odatasets_unique)):
        print(f"\n--- Executing BAT for dataset at position: {i + 1} ---", flush=True)
        try:
            dataset      = odatasets_unique["dataset"].iloc[i]
            dataset_name = odatasets_unique["name"].iloc[i]
            target_raw   = odatasets_unique["class_distribution_vector"].iloc[i]

            target = to_int_list(target_raw)
            if target is None or len(target) < 2:
                print("Target cardinality not available/invalid. Skipping.")
                continue

            row = run_clustering_row(dataset, target, dataset_name)

            if row is not None:
                results.append(row)
                print(f"Time: {row['Execution_Time']:.4f} seconds")
            else:
                print("Dataset skipped.")

        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error processing dataset at position {i + 1}: {e}")

    if results:
        df_out = pd.DataFrame(results)
        out_path = get_predictions_dir() / "pred_BAT.csv"
        df_out.to_csv(out_path, index=False)
        print(f"\nBAT Algorithm finalizado. Archivo guardado en: {out_path}")
    else:
        print("\nBAT Algorithm finalizado sin resultados válidos.")
