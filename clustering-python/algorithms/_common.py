"""
algorithms/_common.py
=====================
Helpers compartidos por todos los algoritmos de clustering.

Funciones públicas:
  - make_numeric_X         : convierte un DataFrame a numérico puro
                             (con SHORTCUT para embeddings: si ya es todo float,
                              salta la lógica de Categorical → 50× más rápido)
  - prepare_data           : separa X (numérico) e y (target categórico)
  - tabulate_clusters      : equivalente vectorizado a tabulate(arr, nbins=k) de R
  - to_int_list            : valida y convierte target_cardinality a list[int]
  - compute_metrics        : ARI, AMI, NMI, Silhouette (coseno)
  - pam_kaufman_rousseeuw  : PAM clásico (BUILD + SWAP), determinístico
  - calculate_centroids    : centroides por media con relleno de clusters vacíos
  - adjust_cardinality     : reasignación al centroide más cercano respetando target
  - evaluate_solution      : silhouette - penalty (función objetivo común)
"""

from __future__ import annotations

import logging
from typing import Any, Optional, cast

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
# sklearn declara silhouette_samples con tipo de retorno ArrayLike (unión muy
# amplia); silenciamos el warning de Pylance.
from sklearn.metrics import (  # pyright: ignore[reportUnknownVariableType]
    adjusted_mutual_info_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_samples,
)

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------------
# Conteo por cluster (equivalente a R tabulate)
# ----------------------------------------------------------------------------

def tabulate_clusters(arr: np.ndarray, k: int) -> np.ndarray:
    """
    Equivalente a `tabulate(arr, nbins=k)` de R: cuenta ocurrencias de los
    valores 1..k en `arr`, ignorando valores fuera de ese rango.
    """
    valid = arr[(arr >= 1) & (arr <= k)]
    if len(valid) == 0:
        return np.zeros(k, dtype=int)
    counts = np.bincount(valid - 1, minlength=k)
    return counts[:k]


# ----------------------------------------------------------------------------
# Preparación de datos
# ----------------------------------------------------------------------------

def _is_pure_float_df(df: pd.DataFrame) -> bool:
    """
    True si todas las columnas son numéricas de tipo float (caso embeddings).
    En ese caso podemos saltar toda la lógica de Categorical y conversión.
    """
    return all(pd.api.types.is_float_dtype(df[c]) for c in df.columns)


def make_numeric_X(x_df: pd.DataFrame) -> Optional[pd.DataFrame]:
    """
    Convierte X a un DataFrame puramente numérico.

    SHORTCUT EMBEDDINGS: si todas las columnas son float (caso embeddings de
    modelos CLIP/BLIP/Qwen-VL/InternVL2), salta toda la lógica de Categorical
    y descarte de columnas constantes (los embeddings no tienen columnas
    constantes en la práctica). Esto es ~50× más rápido para D=512+.

    Para datasets clásicos (mixtos numérico/categórico):
      - Numéricas/booleanas: cast directo a float
      - Factores/categóricas: códigos enteros (1..k)
      - Strings: tratadas como categóricas
      - Descarta columnas no convertibles y columnas constantes
    """
    x_df = pd.DataFrame(x_df)

    # ===== SHORTCUT EMBEDDINGS =====
    if _is_pure_float_df(x_df) and x_df.shape[1] >= 32:
        # Para embeddings reales (D >= 32), confiamos en que no hay columnas
        # constantes ni problemas de tipo. Conversión directa.
        if x_df.shape[1] == 0:
            return None
        return x_df.astype(float, copy=False)

    # ===== RUTA NORMAL (datasets OpenML / tabulares mixtos) =====
    cols: dict[str, np.ndarray] = {}

    for name in x_df.columns:
        col = x_df[name]
        if pd.api.types.is_bool_dtype(col):
            cols[str(name)] = col.astype(float).to_numpy()
        elif pd.api.types.is_numeric_dtype(col):
            cols[str(name)] = col.astype(float).to_numpy()
        elif isinstance(col.dtype, pd.CategoricalDtype):
            cols[str(name)] = col.cat.codes.astype(float).to_numpy() + 1.0
        elif pd.api.types.is_object_dtype(col) or pd.api.types.is_string_dtype(col):
            codes = pd.Categorical(col).codes.astype(float) + 1.0
            cols[str(name)] = codes

    if not cols:
        return None

    out = pd.DataFrame(cols)

    # Descartar columnas constantes (un único valor)
    keep = [c for c in out.columns if out[c].nunique(dropna=True) > 1]
    out = out[keep]
    if out.shape[1] == 0:
        return None
    return out


def prepare_data(dataset: pd.DataFrame) -> Optional[dict[str, Any]]:
    """
    Separa X (features numéricas) e y (target categórico).
    Equivalente a `prepare_data` de R.
    """
    df = pd.DataFrame(dataset)
    if df.shape[1] < 2:
        return None

    if "class" in df.columns:
        y_raw = df["class"]
        x_raw = df.drop(columns=["class"])
    else:
        y_raw = df.iloc[:, -1]
        x_raw = df.iloc[:, :-1]

    x_num = make_numeric_X(x_raw)
    if x_num is None:
        return None

    y_factor = pd.Categorical(y_raw)
    return {"X": x_num, "y": y_factor}


# ----------------------------------------------------------------------------
# Validación de target_cardinality
# ----------------------------------------------------------------------------

def to_int_list(raw: Any) -> Optional[list[int]]:
    """
    Convierte un objeto iterable (lista, ndarray, Series) a list[int].
    Devuelve None si no es convertible o contiene NaN.
    """
    if raw is None:
        return None
    try:
        return [int(c) for c in raw]
    except (TypeError, ValueError):
        return None


# ----------------------------------------------------------------------------
# Métricas
# ----------------------------------------------------------------------------

def compute_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    x_data: np.ndarray,
    dist_matrix: Optional[np.ndarray] = None,
) -> dict[str, float]:
    """
    Calcula ARI, AMI, NMI y Silhouette media (con distancia coseno).
    Cualquier métrica que falle se devuelve como NaN, igual que tryCatch de R.

    Si `dist_matrix` se proporciona (matriz N×N cuadrada), se reutiliza para
    el silhouette en lugar de recalcular pdist.
    """
    out: dict[str, float] = {
        "ARI": float("nan"),
        "AMI": float("nan"),
        "NMI": float("nan"),
        "Silhouette_mean": float("nan"),
    }

    try:
        out["ARI"] = float(adjusted_rand_score(y_true, y_pred))
    except ValueError:
        pass
    try:
        out["AMI"] = float(adjusted_mutual_info_score(y_true, y_pred))
    except ValueError:
        pass
    try:
        out["NMI"] = float(normalized_mutual_info_score(y_true, y_pred))
    except ValueError:
        pass
    try:
        if len(np.unique(y_pred)) >= 2:
            if dist_matrix is None:
                condensed = cast(np.ndarray, pdist(x_data, metric="cosine"))  # pyright: ignore[reportArgumentType, reportCallIssue]
                d_sq = squareform(condensed)
            else:
                d_sq = dist_matrix
            sil_raw = silhouette_samples(d_sq, y_pred, metric="precomputed")
            sil = np.asarray(sil_raw, dtype=float)
            out["Silhouette_mean"] = float(np.nanmean(sil))
    except ValueError:
        pass

    return out


# ----------------------------------------------------------------------------
# PAM (Partitioning Around Medoids) - Kaufman & Rousseeuw (1990)
# ----------------------------------------------------------------------------

def pam_kaufman_rousseeuw(
    d: np.ndarray, k: int, max_iter: int = 100
) -> np.ndarray:
    """
    PAM clásico: BUILD + SWAP. Algoritmo determinístico.

    d : matriz cuadrada de distancias N×N.
    k : número de medoids a encontrar.
    Devuelve un ndarray con los k índices de los medoids (0-based).
    """
    n = d.shape[0]

    # ===== BUILD phase =====
    medoids: list[int] = [int(np.argmin(d.sum(axis=1)))]

    for _ in range(k - 1):
        current_min_d = np.min(d[:, medoids], axis=1)
        best_gain      = -np.inf
        best_candidate = -1
        for i in range(n):
            if i in medoids:
                continue
            gain = float(np.sum(np.maximum(0.0, current_min_d - d[:, i])))
            if gain > best_gain:
                best_gain      = gain
                best_candidate = i
        medoids.append(best_candidate)

    # ===== SWAP phase =====
    for _ in range(max_iter):
        current_cost = float(np.sum(np.min(d[:, medoids], axis=1)))
        best_swap: Optional[tuple[int, int]] = None
        best_cost = current_cost

        for j_idx in range(k):
            for i in range(n):
                if i in medoids:
                    continue
                new_medoids        = medoids.copy()
                new_medoids[j_idx] = i
                new_cost = float(np.sum(np.min(d[:, new_medoids], axis=1)))
                if new_cost < best_cost:
                    best_cost = new_cost
                    best_swap = (j_idx, i)

        if best_swap is None:
            break
        medoids[best_swap[0]] = best_swap[1]

    return np.array(medoids, dtype=int)


# ----------------------------------------------------------------------------
# Centroides y ajuste de cardinalidad (compartidos por BAT, ACO, PSO, etc.)
# ----------------------------------------------------------------------------

def calculate_centroids(
    x_mat: np.ndarray, cluster_assignment: np.ndarray, k: int
) -> np.ndarray:
    """
    Calcula centroides como la media por cluster (etiquetas 1..k).
    Si algún cluster queda vacío, lo rellena con la media global.
    """
    x_mat = np.asarray(x_mat, dtype=float)
    centroids = np.full((k, x_mat.shape[1]), np.nan, dtype=float)

    for c in range(1, k + 1):
        mask = cluster_assignment == c
        if mask.any():
            centroids[c - 1] = x_mat[mask].mean(axis=0)

    if np.isnan(centroids).any():
        global_mu = np.nanmean(x_mat, axis=0)
        nan_rows = np.isnan(centroids).any(axis=1)
        centroids[nan_rows] = global_mu

    return centroids


def adjust_cardinality(
    cluster_assignment: np.ndarray,
    x_mat: np.ndarray,
    centroids: np.ndarray,
    target_cardinality: np.ndarray,
    max_iterations: int = 1000,
) -> np.ndarray:
    """
    Reasigna elementos de clusters sobre-poblados al cluster con espacio
    cuyo centroide queda más cercano (distancia euclídea al cuadrado).
    """
    cluster_assignment = cluster_assignment.copy().astype(int)
    x_mat              = np.asarray(x_mat, dtype=float)
    target             = np.asarray(target_cardinality, dtype=int)
    k                  = len(target)

    cluster_sizes = tabulate_clusters(cluster_assignment, k)
    iteration = 0

    while np.any(cluster_sizes > target) and iteration < max_iterations:
        iteration += 1
        for j in np.where(cluster_sizes > target)[0]:
            cluster_label = int(j) + 1  # nuestras etiquetas son 1-based como en R
            idx = np.where(cluster_assignment == cluster_label)[0]
            if len(idx) == 0:
                continue
            element = idx[0]

            diff      = centroids - x_mat[element]
            distances = np.sum(diff ** 2, axis=1)

            available     = np.where(cluster_sizes < target)[0]
            valid_clusters = available[available != j]

            if len(valid_clusters) == 0:
                cluster_assignment[element] = -1
                cluster_sizes[j] -= 1
            else:
                chosen = valid_clusters[np.argmin(distances[valid_clusters])]
                cluster_assignment[element] = int(chosen) + 1
                cluster_sizes[j]      -= 1
                cluster_sizes[chosen] += 1

    # reasignación final de elementos marcados como -1
    unassigned = np.where(cluster_assignment == -1)[0]
    for element in unassigned:
        available = np.where(cluster_sizes < target)[0]
        if len(available) > 0:
            chosen = available[np.argmin(cluster_sizes[available])]
            cluster_assignment[element] = int(chosen) + 1
            cluster_sizes[chosen] += 1
        else:
            cluster_assignment[element] = 1
            cluster_sizes[0] += 1

    if iteration >= max_iterations:
        logger.warning(
            "adjust_cardinality alcanzó el máximo de iteraciones (%d).",
            max_iterations,
        )

    return cluster_assignment


# ----------------------------------------------------------------------------
# Función objetivo común: silhouette - penalty
# ----------------------------------------------------------------------------

def evaluate_solution(
    cluster_assignment: np.ndarray,
    d_cosine_square: np.ndarray,
    target_cardinality: np.ndarray,
    penalty_weight: float = 10.0,
) -> float:
    """
    Función objetivo: silhouette media (con distancia coseno precomputada)
    menos un penalty proporcional a la desviación de la cardinalidad objetivo.

    Usada por BAT (penalty_weight=10), ACO (penalty_weight=100), PSO, etc.
    """
    if len(np.unique(cluster_assignment)) < 2:
        return float("-inf")

    try:
        sil_raw = silhouette_samples(
            d_cosine_square, cluster_assignment, metric="precomputed"
        )
        sil = np.asarray(sil_raw, dtype=float)
    except ValueError:
        return float("-inf")

    counts  = tabulate_clusters(cluster_assignment, len(target_cardinality))
    penalty = penalty_weight * float(np.sum(np.abs(counts - target_cardinality)))

    # Equivale a mean(ss[, "sil_width"]) - penalty (sin na.rm como en R)
    return float(sil.mean() - penalty)
