"""
embeddings_loader.py
====================
Carga embeddings desde archivos .npz y construye un DataFrame
`odatasets_unique` compatible con la pipeline existente (los 8 algoritmos
y el testing.py).

Estructura esperada de archivos en disco:
    embeddings/
        clip/
            bbc_sport.npz
            ag_news.npz
            ...
        blip/
            ...
        qwen_vl/
            ...
        internvl/
            ...

Cada .npz contiene:
    X            : matriz de embeddings (N, D), float32, idealmente L2-normalizada
    y            : labels enteros (N,)
    class_names  : array de strings con los nombres originales de las clases
    metadata     : dict serializado con info (n, n_classes, source, fecha, etc.)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(message)s")

PROJECT_ROOT          = Path(__file__).resolve().parent
DEFAULT_EMBEDDINGS_DIR = PROJECT_ROOT / "embeddings"

VALID_MODELS = {"clip", "blip", "qwen_vl", "internvl"}


# ============================================================================
# Carga de un único .npz
# ============================================================================

def _load_single_npz(npz_path: Path) -> Optional[dict[str, Any]]:
    """Carga un .npz y devuelve dict con X, y, n, n_classes, distribution."""
    try:
        npz = np.load(npz_path, allow_pickle=True)
    except (OSError, ValueError) as e:
        logger.warning("No se pudo cargar %s: %s", npz_path.name, e)
        return None

    if "X" not in npz or "y" not in npz:
        logger.warning("Archivo %s no tiene claves X/y. Skipping.", npz_path.name)
        return None

    x_emb: np.ndarray = np.asarray(npz["X"], dtype=np.float32)
    y_arr: np.ndarray = np.asarray(npz["y"]).astype(int)

    if x_emb.shape[0] != y_arr.shape[0]:
        logger.warning(
            "Inconsistencia en %s: X tiene %d filas, y tiene %d. Skipping.",
            npz_path.name, x_emb.shape[0], y_arr.shape[0],
        )
        return None

    class_names_raw = npz.get("class_names", None)
    if class_names_raw is None:
        class_names: list[str] = [str(c) for c in sorted(np.unique(y_arr).tolist())]
    else:
        class_names = [str(c) for c in np.asarray(class_names_raw).tolist()]

    unique_classes, counts = np.unique(y_arr, return_counts=True)
    class_distribution_vector = counts.tolist()

    return {
        "dataset_name":              npz_path.stem,
        "X":                         x_emb,
        "y":                         y_arr,
        "class_names":               class_names,
        "n":                         int(x_emb.shape[0]),
        "n_classes":                 int(len(unique_classes)),
        "class_distribution_vector": class_distribution_vector,
    }


# ============================================================================
# Construcción de DataFrame por dataset (Opción A: features + columna class)
# ============================================================================

def _build_dataset_df(x_emb: np.ndarray, y_arr: np.ndarray) -> pd.DataFrame:
    """DataFrame con D columnas de embedding + columna 'class' (string)."""
    d = int(x_emb.shape[1])
    feature_cols = [f"e{i}" for i in range(d)]
    df = pd.DataFrame(x_emb, columns=feature_cols)
    df["class"] = pd.Series(y_arr).astype(str)
    return df


# ============================================================================
# Loader principal
# ============================================================================

def load_all_datasets(
    model: str = "clip",
    embeddings_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Carga todos los embeddings disponibles para un modelo dado."""
    if model not in VALID_MODELS:
        raise ValueError(
            f"Modelo '{model}' no válido. Opciones: {sorted(VALID_MODELS)}"
        )

    base_dir = (embeddings_dir or DEFAULT_EMBEDDINGS_DIR) / model
    if not base_dir.exists():
        raise FileNotFoundError(
            f"Directorio de embeddings no encontrado: {base_dir}\n"
            f"Asegúrate de haber descargado los .npz al directorio correcto."
        )

    npz_files = sorted(base_dir.glob("*.npz"))
    if not npz_files:
        raise FileNotFoundError(
            f"No se encontraron archivos .npz en {base_dir}."
        )

    rows: list[dict[str, Any]] = []
    for npz_path in npz_files:
        info = _load_single_npz(npz_path)
        if info is None:
            continue

        dataset_df = _build_dataset_df(info["X"], info["y"])

        rows.append({
            "name":                      info["dataset_name"],
            "NumberOfClasses":           info["n_classes"],
            "NumberOfInstances":         info["n"],
            "class_distribution_vector": info["class_distribution_vector"],
            "dataset":                   dataset_df,
        })

        logger.info(
            "[%s] %s: N=%d, D=%d, K=%d",
            model, info["dataset_name"],
            info["n"], info["X"].shape[1], info["n_classes"],
        )

    if not rows:
        raise RuntimeError(
            f"Ningún archivo .npz en {base_dir} se pudo cargar correctamente."
        )

    odatasets_unique = pd.DataFrame(rows)
    logger.info(
        "✓ Cargados %d datasets con modelo '%s' desde %s",
        len(odatasets_unique), model, base_dir,
    )
    return odatasets_unique


# ============================================================================
# Carga selectiva
# ============================================================================

def load_one_dataset(
    dataset_name: str,
    model: str = "clip",
    embeddings_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Carga un único dataset por nombre."""
    if model not in VALID_MODELS:
        raise ValueError(
            f"Modelo '{model}' no válido. Opciones: {sorted(VALID_MODELS)}"
        )

    base_dir = (embeddings_dir or DEFAULT_EMBEDDINGS_DIR) / model
    npz_path = base_dir / f"{dataset_name}.npz"

    if not npz_path.exists():
        raise FileNotFoundError(f"No existe {npz_path}")

    info = _load_single_npz(npz_path)
    if info is None:
        raise RuntimeError(f"No se pudo cargar {npz_path}")

    dataset_df = _build_dataset_df(info["X"], info["y"])
    return pd.DataFrame([{
        "name":                      info["dataset_name"],
        "NumberOfClasses":           info["n_classes"],
        "NumberOfInstances":         info["n"],
        "class_distribution_vector": info["class_distribution_vector"],
        "dataset":                   dataset_df,
    }])


# ============================================================================
# CLI para inspección rápida
# ============================================================================

if __name__ == "__main__":
    import sys

    model_arg = sys.argv[1] if len(sys.argv) > 1 else "clip"
    print(f"\nCargando embeddings con modelo '{model_arg}'...")

    try:
        df = load_all_datasets(model=model_arg)
        print("\n=== Resumen ===")
        print(df[["name", "NumberOfClasses", "NumberOfInstances",
                  "class_distribution_vector"]].to_string(index=False))
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}")
        sys.exit(1)
