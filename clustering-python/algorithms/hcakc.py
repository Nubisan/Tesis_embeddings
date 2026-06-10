"""
algorithms/hcakc.py
===================
HCAKC (Hierarchical Clustering Algorithm based on K-means with Constraints).

Implementación fiel al algoritmo del paper:
    Hang, Zhang, Ren & Hu (2009). "A Hierarchical Clustering Algorithm based on
    K-means with Constraints (HCAKC)". 2009 Fourth International Conference on
    Innovative Computing, Information and Control (ICICIC), IEEE, pp. 1479-1482.

------------------------------------------------------------------------------
POR QUÉ ESTA VERSIÓN ES DISTINTA DE LA HEREDADA
------------------------------------------------------------------------------
La versión migrada desde el R original calculaba `IS_values` y una matriz
`Xmat`, pero NUNCA las usaba: la salida final era la simple inicialización por
bloques. Es decir, no ejecutaba el clustering jerárquico aglomerativo que es el
corazón de HCAKC. Esta versión sí lo ejecuta, siguiendo el pseudocódigo de la
Sección III del paper:

    Paso 2     : sobre-segmentación con K-means en t >> K sub-clusters.
    Paso 8     : matriz de cohesión X entre los t sub-clusters (fórmulas 1-3).
    Pasos 9-12 : CUCMC + factor de penalización (mecanismo conservado;
                 con M = C = 0 no modifica X).
    Pasos 13-16: FUSIÓN AGLOMERATIVA — fusiona iterativamente el par de
                 sub-clusters de cohesión máxima hasta dejar K clusters.

------------------------------------------------------------------------------
ADAPTACIÓN A RESTRICCIONES DE CARDINALIDAD (variante de esta tesis)
------------------------------------------------------------------------------
El HCAKC del paper NO impone tamaños de cluster: descubre K con el algoritmo
`Find-K` (máximo de la curva del IS promedio) y sus restricciones son pairwise
(must-link / cannot-link). Esta tesis trabaja con *restricciones de tamaño*
(cardinalidad especificada). Por eso se hacen dos adaptaciones, ambas
documentadas explícitamente:

    1. Se OMITE `Find-K`. El número de clusters K y los tamaños objetivo vienen
       dados por `cluster_sizes` (la cardinalidad objetivo del dataset). El IS
       promedio se sigue calculando y reportando como indicador de calidad,
       pero ya no decide K.
    2. Tras la fusión aglomerativa del paper, se añade un paso final
       `_adjust_cardinality` que reasigna los puntos de frontera para que los K
       clusters cumplan EXACTAMENTE los tamaños objetivo. Es el mismo mecanismo
       usado en los demás algoritmos de la tesis (p. ej. ACO), para mantener
       consistencia metodológica entre algoritmos.

------------------------------------------------------------------------------
NOTA SOBRE LA MÉTRICA DE DISTANCIA
------------------------------------------------------------------------------
El paper define sus fórmulas (1)-(4) con distancia euclidiana. Esta tesis usa
distancia coseno como medida principal (embeddings). El parámetro
`distance_method` controla esto y por defecto es "cosine" para ser coherente con
el resto de los algoritmos de la tesis; cámbialo a "euclidean" si quieres
reproducir literalmente el paper.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

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
# Hiperparámetros
# ============================================================================

@dataclass
class HcakcHyperparams:
    """Hiperparámetros del algoritmo HCAKC (variante con cardinalidad)."""

    algorithm: str = "HCAKC"
    distance_method: str = "cosine"          # "cosine" (tesis) | "euclidean" (paper)
    # t = oversegment_factor * K  (acotado a [K, n]); el paper exige t >> K.
    oversegment_factor: int = 10
    random_state: int = 42
    kmeans_n_init: int = 10
    cohesion_formula: str = "paper_eq_1_2_3_join_exp"
    merge_rule: str = "agglomerative_max_cohesion_group_average"
    constraint_matrices: str = "M_zero_C_zero"
    cardinality_enforcement: str = "adjust_cardinality_post_merge"


HCAKC_HYPERPARAMS: dict[str, Any] = {
    "global": HcakcHyperparams(),
    "runs": {},
}


# ============================================================================
# Bloques geométricos del paper
# ============================================================================

def _dist(a: np.ndarray, b: np.ndarray, metric: str) -> np.ndarray:
    """Distancias fila-a-fila entre a (n×p) y b (m×p) con la métrica dada."""
    return np.asarray(cdist(a, b, metric=metric))  # type: ignore[call-overload]


def _centroids_and_radii(
    data: np.ndarray, labels: np.ndarray, t: int, metric: str
) -> tuple[np.ndarray, np.ndarray]:
    """Centroides y radios de los t sub-clusters.

    El radio sigue la fórmula (1) del paper:
        r(C) = sqrt( (1/|C|) * sum_i d(p_i, c)^2 )
    generalizada a la métrica elegida (euclidiana en el paper, coseno en tesis).
    """
    p = data.shape[1]
    centroids = np.zeros((t, p), dtype=float)
    radii = np.zeros(t, dtype=float)
    for i in range(t):
        members = data[labels == i]
        if members.shape[0] == 0:
            continue
        c = members.mean(axis=0)
        centroids[i] = c
        d = _dist(members, c[None, :], metric).ravel()
        radii[i] = np.sqrt(np.mean(d ** 2))
    return centroids, radii


def _cohesion_matrix(
    data: np.ndarray, labels: np.ndarray, t: int, metric: str
) -> np.ndarray:
    """Matriz de cohesión t x t entre sub-clusters (fórmulas 2 y 3 del paper).

    join(p_i, C_j) = exp( -( d(p_i, c_j) - d(p_i, c_own) ) / r_own )     (eq. 2)
    Chs(C_i, C_j)  = [ sum_{p in C_i} join(p, C_j)
                       + sum_{p in C_j} join(p, C_i) ] / (|C_i| + |C_j|)  (eq. 3)
    """
    n = data.shape[0]
    centroids, radii = _centroids_and_radii(data, labels, t, metric)

    dist = _dist(data, centroids, metric)        # n x t : d(p, c_j)
    d_own = dist[np.arange(n), labels]           # d(p, c_own)
    r_own = radii[labels]
    r_safe = np.where(r_own > 0.0, r_own, 1.0)   # evita división por cero

    expo = np.clip(-(dist - d_own[:, None]) / r_safe[:, None], -50.0, 50.0)
    joinmat = np.exp(expo)                       # n x t

    s_mat = np.zeros((t, t), dtype=float)
    sizes = np.zeros(t, dtype=float)
    for i in range(t):
        members = labels == i
        sizes[i] = members.sum()
        if sizes[i] > 0:
            s_mat[i, :] = joinmat[members, :].sum(axis=0)

    denom = sizes[:, None] + sizes[None, :]
    denom = np.where(denom > 0.0, denom, 1.0)
    return (s_mat + s_mat.T) / denom             # eq. 3


def _cucmc(
    chs: np.ndarray,
    must_link: list[tuple[int, int]],
    cannot_link: list[tuple[int, int]],
) -> np.ndarray:
    """CUCMC: Constraint-based Update of Cohesion Matrix between Clusters.

    Definición 2 del paper. Con M = C = 0 (caso de esta tesis) no modifica nada;
    el mecanismo se conserva por fidelidad y para soportar restricciones futuras.
    """
    chs = chs.copy()
    for i, j in must_link:
        neigh = np.maximum(chs[i, :], chs[j, :])   # propagación a vecinos
        chs[i, :] = neigh
        chs[:, i] = neigh
        chs[j, :] = neigh
        chs[:, j] = neigh
        chs[i, j] = chs[j, i] = 1.0
    for i, j in cannot_link:
        chs[i, j] = chs[j, i] = 0.0
    return chs


def _agglomerative_merge(
    data: np.ndarray,
    labels: np.ndarray,
    k_target: int,
    metric: str,
    must_link: list[tuple[int, int]],
    cannot_link: list[tuple[int, int]],
) -> np.ndarray:
    """Pasos 8-16 del paper: fusiona el par de cohesión máxima hasta K clusters.

    La matriz de cohesión X se calcula UNA sola vez sobre los t micro-clusters
    (paso 8). Luego, en cada paso se extrae el par (C_i, C_j) de cohesión máxima,
    se fusionan y t := t - 1, hasta que t = K (pasos 13-16).

    Al fusionar C_i y C_j, la cohesión del cluster resultante con cualquier otro
    C_r se recombina con el PROMEDIO PONDERADO POR TAMAÑO:
        Chs(C_i|C_j, C_r) := (|C_i|·Chs(C_i,C_r) + |C_j|·Chs(C_j,C_r)) / (|C_i|+|C_j|)
    Esta actualización es consistente con la fórmula (3) del paper (que ya es un
    promedio de `join` normalizado por tamaño): equivale a recalcular la cohesión
    del cluster unión manteniendo fijos los `join` de los micro-clusters
    originales. Es el enlace "group-average" (UPGMA); evita tanto el efecto bola
    de nieve de recalcular centroides como el encadenamiento del enlace simple.
    """
    labels = labels.copy()
    active = sorted(set(labels.tolist()))
    remap = {old: new for new, old in enumerate(active)}
    labels = np.array([remap[v] for v in labels], dtype=int)
    t = len(active)
    if t <= k_target:
        return labels

    chs = _cohesion_matrix(data, labels, t, metric)   # paso 8 (una vez)
    chs = _cucmc(chs, must_link, cannot_link)         # no-op con M = C = 0

    sim = chs.copy()
    np.fill_diagonal(sim, -np.inf)                    # nunca fusionar consigo mismo

    members: dict[int, list[int]] = {i: [i] for i in range(t)}
    csize = np.array([int(np.sum(labels == i)) for i in range(t)], dtype=float)
    alive: list[int] = list(range(t))

    while len(alive) > k_target:                      # pasos 13-16
        sub = np.array(alive)
        sub_sim = sim[np.ix_(sub, sub)]
        ai, aj = np.unravel_index(int(np.argmax(sub_sim)), sub_sim.shape)
        i, j = int(sub[ai]), int(sub[aj])
        if i == j or not np.isfinite(sim[i, j]):
            break

        members[i].extend(members[j])
        wi, wj = csize[i], csize[j]
        new_row = (wi * sim[i, :] + wj * sim[j, :]) / (wi + wj)
        sim[i, :] = new_row
        sim[:, i] = new_row
        sim[i, i] = -np.inf
        csize[i] = wi + wj
        sim[j, :] = -np.inf
        sim[:, j] = -np.inf
        alive.remove(j)

    micro_to_super: dict[int, int] = {}
    for super_id in alive:
        for micro_id in members[super_id]:
            micro_to_super[micro_id] = super_id
    new_labels = np.array([micro_to_super[v] for v in labels], dtype=int)

    remap2 = {old: new for new, old in enumerate(sorted(set(new_labels.tolist())))}
    return np.array([remap2[v] for v in new_labels], dtype=int)


def _ensure_k_clusters(
    labels: np.ndarray, data: np.ndarray, k: int, random_state: int
) -> np.ndarray:
    """Salvaguarda: si la fusión dejó < K clusters, parte el mayor con 2-means."""
    labels = labels.copy()
    uniq = sorted(set(labels.tolist()))
    while len(uniq) < k:
        sizes = {c: int(np.sum(labels == c)) for c in uniq}
        biggest = max(sizes, key=sizes.get)
        idx = np.where(labels == biggest)[0]
        if len(idx) < 2:
            break
        from sklearn.cluster import KMeans
        sub = KMeans(n_clusters=2, n_init=10, random_state=random_state).fit_predict(data[idx])
        labels[idx[sub == 1]] = max(uniq) + 1
        uniq = sorted(set(labels.tolist()))
    remap = {old: new for new, old in enumerate(sorted(set(labels.tolist())))}
    return np.array([remap[v] for v in labels], dtype=int)


def _is_vectorized(
    data: np.ndarray, clusters: np.ndarray, centroids: np.ndarray, metric: str
) -> np.ndarray:
    """IS (Improved Silhouette) del paper, fórmula (4).

        IS(o_i) = (b_i - a_i) / max(a_i, b_i)
    a_i = distancia al centroide propio; b_i = mínima distancia a los demás
    centroides. `clusters` es 1-based (como en el resto de la tesis).
    """
    n = data.shape[0]
    cluster_idx = (clusters - 1).astype(int)
    d_mat = _dist(data, centroids, metric)
    a_arr = d_mat[np.arange(n), cluster_idx]
    d_other = d_mat.copy()
    d_other[np.arange(n), cluster_idx] = np.inf
    b_arr = d_other.min(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        is_values = (b_arr - a_arr) / np.maximum(a_arr, b_arr)
    return np.nan_to_num(is_values, nan=0.0)


# ============================================================================
# Restricción de cardinalidad (mismo mecanismo que ACO)
# ============================================================================

def _adjust_cardinality(
    labels: np.ndarray, data: np.ndarray, target_sizes: np.ndarray, metric: str
) -> np.ndarray:
    """Reasigna puntos de frontera hasta que cada cluster cumpla su tamaño.

    Las etiquetas resultantes de la fusión son arbitrarias, así que primero se
    emparejan los tamaños objetivo con los clusters por rango (el objetivo más
    grande al cluster más grande, etc.) para minimizar los puntos que se mueven.
    `labels` aquí es 0-based.
    """
    labels = labels.copy()
    k = len(target_sizes)
    target = np.asarray(target_sizes, dtype=int)

    cur_sizes = np.array([int(np.sum(labels == c)) for c in range(k)])
    order_clusters = np.argsort(cur_sizes)[::-1]
    order_targets = np.argsort(target)[::-1]
    target_for_cluster = np.zeros(k, dtype=int)
    for rank in range(k):
        target_for_cluster[order_clusters[rank]] = target[order_targets[rank]]

    sizes = cur_sizes.copy()
    max_iter, it = max(1000, 2 * len(labels)), 0
    while np.any(sizes > target_for_cluster) and it < max_iter:
        it += 1
        centroids, _ = _centroids_and_radii(data, labels, k, metric)
        for j in np.where(sizes > target_for_cluster)[0]:
            idx = np.where(labels == j)[0]
            if len(idx) == 0:
                continue
            element = idx[0]
            dists = _dist(centroids, data[element][None, :], metric).ravel()
            available = np.where(sizes < target_for_cluster)[0]
            valid = available[available != j]
            if len(valid) == 0:
                labels[element] = -1
                sizes[j] -= 1
            else:
                chosen = valid[int(np.argmin(dists[valid]))]
                labels[element] = chosen
                sizes[j] -= 1
                sizes[chosen] += 1

    for element in np.where(labels == -1)[0]:
        available = np.where(sizes < target_for_cluster)[0]
        chosen = available[int(np.argmin(sizes[available]))] if len(available) else 0
        labels[element] = chosen
        sizes[chosen] += 1

    return labels


# ============================================================================
# HCAKC principal (variante de cardinalidad especificada)
# ============================================================================

def hcakc_specified(
    data: np.ndarray,
    k: int,
    cluster_sizes: np.ndarray,
    distance_method: str = "cosine",
    oversegment_factor: int = 10,
    random_state: int = 42,
    kmeans_n_init: int = 10,
    must_link: Optional[list[tuple[int, int]]] = None,
    cannot_link: Optional[list[tuple[int, int]]] = None,
) -> np.ndarray:
    """Ejecuta HCAKC y devuelve las etiquetas finales (1-based, longitud n).

    Flujo: K-means a t >> K sub-clusters -> cohesión (eq. 1-3) -> CUCMC ->
    fusión aglomerativa hasta K -> ajuste de cardinalidad.
    """
    data = np.asarray(data, dtype=float)
    n = data.shape[0]
    cluster_sizes = np.asarray(cluster_sizes, dtype=int)

    if k < 2:
        raise ValueError("k < 2")
    if int(cluster_sizes.sum()) != n:
        raise ValueError("sum(cluster_sizes) != n")

    must_link = must_link or []
    cannot_link = cannot_link or []

    # ----- Paso 2 del paper: sobre-segmentación con K-means (t >> K) -----
    t0 = min(n, max(2 * k, oversegment_factor * k))
    t0 = max(t0, k)
    from sklearn.cluster import KMeans

    micro = KMeans(
        n_clusters=t0, n_init=kmeans_n_init, random_state=random_state
    ).fit_predict(data)

    active = sorted(set(micro.tolist()))
    remap = {old: new for new, old in enumerate(active)}
    micro = np.array([remap[v] for v in micro], dtype=int)
    t0 = len(active)

    # ----- Pasos 8-16: cohesión + fusión aglomerativa hasta K -----
    if t0 > k:
        merged = _agglomerative_merge(
            data, micro, k, distance_method, must_link, cannot_link
        )
    else:
        merged = micro
    merged = _ensure_k_clusters(merged, data, k, random_state)

    # ----- Adaptación de cardinalidad: tamaños exactos -----
    final = _adjust_cardinality(merged, data, cluster_sizes, distance_method)

    return final + 1  # 1-based, como el resto de la tesis


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

    hp: HcakcHyperparams = HCAKC_HYPERPARAMS["global"]
    return hcakc_specified(
        x_mat,
        k,
        target,
        distance_method=hp.distance_method,
        oversegment_factor=hp.oversegment_factor,
        random_state=hp.random_state,
        kmeans_n_init=hp.kmeans_n_init,
    )


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

    x_df: pd.DataFrame = data["X"]
    y_factor: pd.Categorical = data["y"]

    target = np.asarray(target_cardinality, dtype=int)
    x_mat = x_df.to_numpy(dtype=float)
    n = int(x_mat.shape[0])

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
    total_time = time.perf_counter() - start_total

    y_int = (y_factor.codes + 1).astype(int)  # 1-based como as.integer(factor) de R
    metrics = compute_metrics(y_int, y_pred, x_mat)

    hp: HcakcHyperparams = HCAKC_HYPERPARAMS["global"]
    centroids, _ = _centroids_and_radii(x_mat, y_pred - 1, len(target), hp.distance_method)
    is_mean = float(np.mean(_is_vectorized(x_mat, y_pred, centroids, hp.distance_method)))

    HCAKC_HYPERPARAMS["runs"][dataset_name] = {
        "k":                       len(target),
        "n":                       len(y_pred),
        "n_features":              x_mat.shape[1],
        "target_sizes":            target.tolist(),
        "distance_method":         hp.distance_method,
        "oversegment_factor":      hp.oversegment_factor,
        "merge_rule":              hp.merge_rule,
        "constraint_matrices":     hp.constraint_matrices,
        "cardinality_enforcement": hp.cardinality_enforcement,
        "IS_mean":                 is_mean,
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
        out_path = get_predictions_dir() / "pred_HCAKC.csv"
        df_out.to_csv(out_path, index=False)
        print(f"\nHCAKC finalizado. Archivo guardado en: {out_path}")
    else:
        print("\nHCAKC finalizado sin resultados válidos. No se generó archivo CSV.")