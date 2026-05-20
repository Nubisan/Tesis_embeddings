"""
algorithms/hcakc.py
===================
HCAKC (Hard Clustering with Algorithmic Cardinality Constraints).
Migración Python del archivo HCAKC.R.

Estructura del algoritmo (con M=0, C=0 — sin must-link / cannot-link):
  1. Inicialización por bloques secuenciales: los primeros target[1] puntos
     van al cluster 1, los siguientes target[2] al cluster 2, etc.
  2. Calcula los centroides como medias de cada bloque.
  3. Calcula el Index Score (IS) para cada punto, similar al Silhouette pero
     usando distancias euclidianas a centroides en lugar de a puntos.
  4. Calcula la matriz Xmat K×K (productos cruzados de IS por cluster).
  5. Aplica CUCMC (con M=0, C=0 no modifica nada).

Nota importante (replicada fielmente desde R):
  El R calcula `IS_values` y `Xmat` pero NUNCA los usa para reasignar puntos.
  La salida final son los `clusters` de la inicialización por bloques. Si los
  datos vienen ordenados por clase (como iris), HCAKC funciona muy bien.
  Si vienen mezclados, da resultados ~aleatorios.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, cast

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from ._common import (
    compute_metrics,
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
class HcakcHyperparams:
    """Hiperparámetros globales del algoritmo HCAKC."""
    algorithm:           str = "HCAKC"
    distance_method:     str = "cosine"
    silhouette_metric:   str = "euclidean"
    constraint_matrices: str = "M_zero_C_zero"
    centroid_update:     str = "colMeans"
    init_method:         str = "sequential_block"
    objective:           str = "IS_cross_product_with_CUCMC"


HCAKC_HYPERPARAMS: dict[str, Any] = {
    "global": HcakcHyperparams(),
    "runs":   {},
}

PROJECT_ROOT    = Path(__file__).resolve().parent.parent
PREDICTIONS_DIR = PROJECT_ROOT / "predictions"
PREDICTIONS_DIR.mkdir(parents=True, exist_ok=True)
PRED_HCAKC_CSV  = PREDICTIONS_DIR / "pred_HCAKC.csv"


# ============================================================================
# Funciones auxiliares HCAKC
# ============================================================================

def _is_vectorized(
    data: np.ndarray, clusters: np.ndarray, centroids: np.ndarray
) -> np.ndarray:
    """
    Calcula el Index Score (IS) para cada punto:
        IS_i = (b_i - a_i) / max(a_i, b_i)
    donde:
        a_i = distancia euclidiana al centroide del propio cluster
        b_i = mínima distancia euclidiana a un centroide de otro cluster

    OPT vs R: usa scipy.spatial.distance.cdist (BLAS) en una sola llamada,
    en lugar del loop manual sobre centroides del R original.
    """
    n = int(data.shape[0])

    # Distancias euclidianas (n × K) en una sola operación BLAS.
    # Bug del stub de scipy: cdist sí acepta strings como `metric`.
    d_mat = cast(np.ndarray, cdist(data, centroids, metric="euclidean"))
    # pyright: ignore[reportArgumentType, reportCallIssue]

    # clusters es 1-based; convertimos a 0-based para indexar la matriz
    cluster_idx = (clusters - 1).astype(int)

    # a_i = distancia al propio centroide
    a_arr = d_mat[np.arange(n), cluster_idx]

    # b_i = mínima distancia a cualquier OTRO centroide
    # Ponemos inf en la columna del propio cluster para excluirlo del min
    d_other = d_mat.copy()
    d_other[np.arange(n), cluster_idx] = np.inf
    b_arr = d_other.min(axis=1)

    # (b - a) / max(a, b). Silenciamos warnings por si max=0 (raro pero posible).
    with np.errstate(divide="ignore", invalid="ignore"):
        is_values = (b_arr - a_arr) / np.maximum(a_arr, b_arr)
    return is_values


def _cucmc(
    xm: np.ndarray, m_mat: np.ndarray, c_mat: np.ndarray
) -> np.ndarray:
    """
    CUCMC (Constraint Update on Cardinality Constraint Matrix).
    Aplica las matrices M (must-link) y C (cannot-link).
    Con M=0 y C=0 (caso de este algoritmo), no modifica nada.
    """
    xm = xm.copy()
    xm[m_mat == 1] = 1.0
    xm[(c_mat == 1) & (m_mat != 1)] = 0.0
    return xm


# ============================================================================
# HCAKC con M=0, C=0
# ============================================================================

def hcakc_specified(
    data: np.ndarray, k: int, cluster_sizes: np.ndarray
) -> np.ndarray:
    """
    HCAKC con M=0, C=0 (sin restricciones must-link / cannot-link).

    NOTA: el R calcula IS_values y Xmat pero NO los usa para reasignar puntos.
    La asignación final es la inicialización secuencial por bloques. Replicamos
    este comportamiento con fidelidad para que los resultados coincidan con R.
    """
    n, p = int(data.shape[0]), int(data.shape[1])

    # ----- Inicialización por bloques secuenciales -----
    clusters: np.ndarray  = np.zeros(n, dtype=int)
    centroids: np.ndarray = np.zeros((k, p), dtype=float)

    start_idx = 0
    for cluster_label in range(1, k + 1):
        end_idx = start_idx + int(cluster_sizes[cluster_label - 1])
        if end_idx > n:
            raise ValueError("Índice fuera de rango en inicialización.")
        clusters[start_idx:end_idx]      = cluster_label  # 1-based
        centroids[cluster_label - 1, :]  = data[start_idx:end_idx, :].mean(axis=0)
        start_idx = end_idx

    # ----- Cálculo de IS values (vectorizado) -----
    # Replicado por fidelidad con R, aunque el resultado no se use para
    # reasignación final. Útil si en el futuro M o C son distintas de cero.
    is_values = _is_vectorized(data, clusters, centroids)

    # ----- Cálculo de Xmat (matriz K×K de productos cruzados) -----
    # También replicado por fidelidad. M y C son ceros con este algoritmo.
    m_mat = np.zeros((k, k), dtype=float)
    c_mat = np.zeros((k, k), dtype=float)

    xmat = np.zeros((k, k), dtype=float)
    for i in range(1, k + 1):
        is_i = is_values[clusters == i]
        if len(is_i) == 0:
            continue
        for j in range(1, k + 1):
            is_j = is_values[clusters == j]
            if len(is_j) == 0:
                continue
            min_len = int(min(len(is_i), len(is_j)))
            xmat[i - 1, j - 1] = float(np.nansum(is_i[:min_len] * is_j[:min_len]))

    xmat = _cucmc(xmat, m_mat, c_mat)
    # xmat se descarta intencionalmente — el algoritmo en R hace lo mismo.
    _ = xmat

    return clusters


def run_hcakc_algo(
    x_mat: np.ndarray,
    y_factor: pd.Categorical,
    target_cardinality: list[int],
) -> np.ndarray:
    """Wrapper con validaciones equivalentes al R original."""
    if len(x_mat) != len(y_factor):
        raise ValueError("nrow(X) != length(y)")
    target = np.asarray(target_cardinality, dtype=int)
    if np.any(np.isnan(target.astype(float))):
        raise ValueError("target_cardinality contiene NA")
    if int(target.sum()) != x_mat.shape[0]:
        raise ValueError("Suma de cardinalidades != nrow(X)")
    k = int(len(target))
    if k < 2:
        raise ValueError("K < 2")

    return hcakc_specified(x_mat, k, target)


# ============================================================================
# Runner por dataset
# ============================================================================

def run_clustering_row(
    dataset: pd.DataFrame,
    target_cardinality: list[int],
    dataset_name: str,
) -> Optional[dict[str, Any]]:
    """Ejecuta HCAKC sobre un dataset y devuelve la fila de resultados."""
    data = prepare_data(dataset)
    if data is None:
        return None

    x_df:     pd.DataFrame    = data["X"]
    y_factor: pd.Categorical = data["y"]

    target = np.asarray(target_cardinality, dtype=int)
    x_mat  = x_df.to_numpy(dtype=float)
    n      = int(x_mat.shape[0])

    if int(target.sum()) != n:
        return None
    if len(target) < 2:
        return None

    start_total = time.perf_counter()

    try:
        y_pred = run_hcakc_algo(x_mat, y_factor, target_cardinality)
    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.warning("Error en run_hcakc_algo para %s: %s", dataset_name, e)
        return None

    y_int      = (y_factor.codes + 1).astype(int)  # 1-based como as.integer(factor) de R
    total_time = time.perf_counter() - start_total

    metrics = compute_metrics(y_int, y_pred, x_mat)

    HCAKC_HYPERPARAMS["runs"][dataset_name] = {
        "k":                   len(target),
        "n":                   len(y_pred),
        "n_features":          x_mat.shape[1],
        "target_sizes":        target.tolist(),
        "distance_method":     "cosine",
        "silhouette_metric":   "euclidean",
        "constraint_matrices": "M_zero_C_zero",
        "init_method":         "sequential_block",
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
    """Loop principal equivalente al de HCAKC.R."""
    results: list[dict[str, Any]] = []

    for i in range(len(odatasets_unique)):
        print(
            f"\n--- Executing HCAKC for dataset at position: {i + 1} ---",
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
            print(f"Error processing dataset at position {i + 1}: {e}")

    if results:
        df_out = pd.DataFrame(results)
        df_out.to_csv(PRED_HCAKC_CSV, index=False)
        print(f"\nHCAKC finalizado. Archivo guardado en: {PRED_HCAKC_CSV}")
    else:
        print("\nHCAKC finalizado sin resultados válidos. No se generó archivo CSV.")
