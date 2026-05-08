"""
algorithms/kmedoids.py
======================
K-MedoidsSC: K-Medoids con Restricciones de Tamaño (Size-Constrained).
Migración Python del archivo Kmedoids.R.

Algoritmo en dos etapas:
  1. PAM (k-medoids) para encontrar k medoids representativos
  2. Asignación greedy ordenada: cada punto se ordena por su distancia al
     medoid más cercano y se asignan los target[i] más cercanos al cluster i,
     en orden 1..k. Garantiza la cardinalidad exacta sin necesidad de ILP.

Diferencias con R:
  - PAM implementado a mano (compartido con CSCLP vía _common.py)
  - Operaciones vectorizadas con NumPy en lugar de bucles R
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Optional, cast

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from ._common import (
    compute_metrics,
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

KMEDOIDS_HYPERPARAMS_GLOBAL: dict[str, Any] = {
    "algorithm":         "K-MedoidsSC",
    "distance_method":   "cosine",
    "medoid_init":       "pam",
    "assignment_method": "sorted_greedy",
}

# Hiperparámetros guardados por dataset durante la ejecución
KMEDOIDS_RUNS: dict[str, dict[str, Any]] = {}

PROJECT_ROOT      = Path(__file__).resolve().parent.parent
PREDICTIONS_DIR   = PROJECT_ROOT / "predictions"
PREDICTIONS_DIR.mkdir(parents=True, exist_ok=True)
PRED_KMEDOIDS_CSV = PREDICTIONS_DIR / "pred_KMEDOIDS.csv"


# ============================================================================
# Asignación con restricción de cardinalidad (sorted greedy)
# ============================================================================

def sc_medoids(
    d: np.ndarray,
    k: int,
    target_sizes: np.ndarray,
    medoid_idx: np.ndarray,
) -> dict[str, Any]:
    """
    Asigna puntos a medoids con cardinalidad fija usando una estrategia greedy:

    1. Para cada punto, calcula su distancia al medoid más cercano (vector min).
    2. Ordena los puntos por esa distancia (ascendente; los más "fáciles" primero).
    3. Asigna los target[i] primeros al cluster i, en orden 1..k.
    4. Si quedaran puntos sobrantes, los asigna a su medoid más cercano.

    Esta es la lógica de SC_medoids en R, optimizada con NumPy.
    """
    n = int(d.shape[0])

    # OPT: D_C = D[:, medoid_idx] una sola vez (en R se calculaba 2 veces)
    d_c = d[:, medoid_idx]  # shape (n, k)

    # Asignación inicial: cada punto a su medoid más cercano (1-based como R)
    cl: np.ndarray = np.argmin(d_c, axis=1).astype(int) + 1

    # Ordenar puntos por su distancia mínima al conjunto de medoids (ascendente).
    # Equivalente a `order(apply(D_C, 1, min))` de R.
    sorted_points = np.argsort(d_c.min(axis=1))

    # Asignación greedy por cluster
    target = np.asarray(target_sizes, dtype=int)
    pos = 0
    for i in range(1, k + 1):
        size_i = int(target[i - 1])
        cl[sorted_points[pos:pos + size_i]] = i
        pos += size_i

    # Si queda algún punto sobrante (sólo si sum(target) < n), lo asignamos al
    # medoid más cercano. En R esto está por seguridad; con sum(target)==n
    # ya validado, este bloque no se ejecuta.
    if pos < n:
        for point in sorted_points[pos:]:
            cl[point] = int(np.argmin(d[point, medoid_idx])) + 1

    return {
        "medoids":    medoid_idx,
        "clustering": cl,
        "y_predict":  cl.astype(int),
    }


# ============================================================================
# Runner por dataset
# ============================================================================

def run_kmedoids_sc_row(
    dataset: pd.DataFrame,
    target_cardinality: list[int],
    dataset_name: str,
) -> Optional[dict[str, Any]]:
    """Ejecuta K-MedoidsSC sobre un dataset y devuelve la fila de resultados."""
    prepared = prepare_data(dataset)
    if prepared is None:
        return None

    x_df:     pd.DataFrame    = prepared["X"]
    y_factor: pd.Categorical = prepared["y"]

    target = np.asarray(target_cardinality, dtype=int)
    x_mat  = x_df.to_numpy(dtype=float)
    n      = int(x_mat.shape[0])
    k      = int(len(target))

    # Validaciones equivalentes a las del R original
    if k < 2:
        return None
    if int(target.sum()) != n:
        return None

    start_total = time.perf_counter()

    # ----- Distancia coseno N×N (reusada para silhouette) -----
    # Bug del stub de scipy: pdist sí acepta strings como `metric`.
    condensed = cast(np.ndarray, pdist(x_mat, metric="cosine"))  # pyright: ignore[reportArgumentType, reportCallIssue]
    full_dist = squareform(condensed)

    # ----- PAM para medoids iniciales -----
    medoid_idx = pam_kaufman_rousseeuw(full_dist, k)

    # ----- Asignación con restricción de cardinalidad -----
    result = sc_medoids(full_dist, k, target, medoid_idx)

    total_time = time.perf_counter() - start_total

    y_pred = np.asarray(result["y_predict"], dtype=int)
    y_int  = (y_factor.codes + 1).astype(int)  # 1-based como as.integer(factor) de R

    metrics = compute_metrics(y_int, y_pred, x_mat, dist_matrix=full_dist)

    # Guardar hiperparámetros de la corrida
    KMEDOIDS_RUNS[dataset_name] = {
        "k":                 k,
        "n":                 n,
        "n_features":        x_mat.shape[1],
        "target_sizes":      target.tolist(),
        "medoid_indices":    medoid_idx.tolist(),
        "distance_method":   "cosine",
        "medoid_init":       "pam",
        "assignment_method": "sorted_greedy",
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
    """Loop principal equivalente al de Kmedoids.R."""
    results: list[dict[str, Any]] = []

    for i in range(len(odatasets_unique)):
        print(
            f"\n--- Executing KMEDOIDS for dataset at position: {i + 1} ---",
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

            row = run_kmedoids_sc_row(dataset, target, dataset_name)

            if row is not None:
                results.append(row)
                print(f"Time: {row['Execution_Time']:.4f} seconds")
            else:
                print("Dataset skipped.")

        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error processing dataset at position {i + 1}: {e}")

    if results:
        df_out = pd.DataFrame(results)
        df_out.to_csv(PRED_KMEDOIDS_CSV, index=False)
        print(
            f"\nK-Medoids finalizado. Archivo guardado en: {PRED_KMEDOIDS_CSV}"
        )
    else:
        print("\nK-Medoids finalizado sin resultados válidos.")
