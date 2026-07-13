"""
embeddings_loader.py
====================
Carga embeddings desde archivos .npz y construye el DataFrame
`odatasets_unique` compatible con la pipeline existente (los 8 algoritmos
y testing.py). Reemplaza a Openml.R / openml_loader.

Estructura esperada en disco (NUEVA: con sub-carpeta de modalidad):
    embeddings/
        <modelo>/                      modelo in {clip, blip, qwen_vl, internvl}
            imagen/
                orl.npz
                mnist.npz
                ...
            texto/
                bbc_sport.npz
                ...
    (retrocompatible: si no existe la sub-carpeta de modalidad, lee los .npz
     directamente de embeddings/<modelo>/ como en la version anterior.)

Cada .npz contiene:
    X            : (N, D) float32, idealmente L2-normalizado
    y            : (N,) labels enteros
    class_names  : (k,) nombres de clase
    grupo        : 'A' | 'B'   (opcional)
    k            : int         (opcional)
    modalidad    : 'imagen' | 'texto'  (opcional)
    modelo       : str         (opcional)

Funcion publica principal (firma compatible con testing.py):
    load_all_datasets(model, modalidad='imagen', grupo=None) -> pd.DataFrame
        Devuelve `odatasets_unique` con columnas:
          - name                       (str)
          - NumberOfClasses            (int)
          - NumberOfInstances          (int)
          - class_distribution_vector  (list[int]  -> es el target_cardinality)
          - dataset                    (pd.DataFrame: D columnas float + columna 'class')
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(message)s")

PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_EMBEDDINGS_DIR = Path(os.environ.get("EMB_ROOT", str(PROJECT_ROOT / "embeddings")))

VALID_MODELS = {"clip", "blip", "qwen_vl", "internvl"}
VALID_MODALIDADES = {"imagen", "texto"}


# ============================================================================
# Carga de un unico .npz -> dict con metadatos + DataFrame interno
# ============================================================================

def _load_single_npz(npz_path: Path) -> Optional[dict[str, Any]]:
    """
    Carga un .npz y devuelve un dict con:
        name, NumberOfClasses, NumberOfInstances,
        class_distribution_vector (list[int]), dataset (DataFrame D float + 'class'),
        grupo.
    Devuelve None si el archivo es invalido (se omite, no rompe la corrida).
    """
    try:
        npz = np.load(npz_path, allow_pickle=True)
    except (OSError, ValueError) as e:
        logger.warning("No se pudo cargar %s: %s", npz_path.name, e)
        return None

    if "X" not in npz or "y" not in npz:
        logger.warning("Archivo %s sin claves X/y. Omitido.", npz_path.name)
        return None

    x_emb = np.asarray(npz["X"], dtype=np.float32)
    y_raw = np.asarray(npz["y"]).astype(int)

    if x_emb.ndim != 2 or x_emb.shape[0] != y_raw.shape[0] or x_emb.shape[0] == 0:
        logger.warning(
            "Inconsistencia en %s: X=%s, y=%s. Omitido.",
            npz_path.name, getattr(x_emb, "shape", None), getattr(y_raw, "shape", None),
        )
        return None

    # Remapear etiquetas a 0..k-1 CONSECUTIVAS (evita el bug 'index out of bounds'
    # cuando un dataset trae labels originales no consecutivos tras la limpieza).
    _, y = np.unique(y_raw, return_inverse=True)
    y = y.astype(int)

    n, _d = x_emb.shape
    k = int(y.max()) + 1
    counts = np.bincount(y, minlength=k).astype(int)  # conteo por clase (= target_cardinality)

    # DataFrame interno: D columnas float + columna 'class' (lo que espera prepare_data)
    inner = pd.DataFrame(x_emb)                       # columnas 0..D-1, dtype float32
    inner.columns = [str(c) for c in inner.columns]   # nombres de columna como str
    inner["class"] = y                                # columna target

    grupo = str(npz["grupo"]) if "grupo" in npz else ""

    return {
        "name": npz_path.stem,
        "NumberOfClasses": int(k),
        "NumberOfInstances": int(n),
        "class_distribution_vector": [int(c) for c in counts],
        "dataset": inner,
        "grupo": grupo,
    }


# ============================================================================
# Resolucion de carpeta (con modalidad y retrocompatibilidad)
# ============================================================================

def _resolver_dir(model: str, modalidad: str, embeddings_dir: Path) -> Path:
    """
    Devuelve embeddings/<model>/<modalidad>/ si existe; si no, cae a
    embeddings/<model>/ (estructura plana antigua).
    """
    base = embeddings_dir / model
    con_modalidad = base / modalidad
    if con_modalidad.is_dir():
        return con_modalidad
    return base


# ============================================================================
# Punto de entrada principal (firma usada por testing.py)
# ============================================================================

def load_all_datasets(
    model: str = "clip",
    modalidad: str = "imagen",
    grupo: Optional[str] = None,
    embeddings_dir: Optional[Path | str] = None,
) -> pd.DataFrame:
    """
    Construye `odatasets_unique` (DataFrame) leyendo todos los .npz de
    embeddings/<model>/<modalidad>/.

    Parameters
    ----------
    model       : 'clip' | 'blip' | 'qwen_vl' | 'internvl'
    modalidad   : 'imagen' | 'texto'      (NUEVO; antes solo existia imagen)
    grupo       : 'A' | 'B' | None        (None = ambos grupos)
    embeddings_dir : raiz de embeddings (default: ./embeddings o $EMB_ROOT)

    Returns
    -------
    pd.DataFrame con columnas:
        name, NumberOfClasses, NumberOfInstances, class_distribution_vector, dataset
    (mismo formato que el loader original; los algoritmos no cambian.)
    """
    if model not in VALID_MODELS:
        raise ValueError(f"model debe ser uno de {sorted(VALID_MODELS)}, recibido: {model!r}")
    if modalidad not in VALID_MODALIDADES:
        raise ValueError(f"modalidad debe ser una de {sorted(VALID_MODALIDADES)}, recibido: {modalidad!r}")

    emb_dir = Path(embeddings_dir) if embeddings_dir is not None else DEFAULT_EMBEDDINGS_DIR
    carpeta = _resolver_dir(model, modalidad, emb_dir)

    if not carpeta.is_dir():
        raise FileNotFoundError(
            f"No existe la carpeta {carpeta}. ¿Descomprimiste emb_{model}.zip "
            f"bajo embeddings/{model}/{modalidad}/?"
        )

    filas: list[dict[str, Any]] = []
    for npz_path in sorted(carpeta.glob("*.npz")):
        info = _load_single_npz(npz_path)
        if info is None:
            continue
        if grupo is not None and info["grupo"] != str(grupo):
            continue
        filas.append(info)

    if not filas:
        raise FileNotFoundError(
            f"No se cargo ningun .npz valido en {carpeta} "
            f"(modelo={model}, modalidad={modalidad}, grupo={grupo})."
        )

    # Mismas 5 columnas y orden que el loader original (los algoritmos las usan por nombre)
    odatasets_unique = pd.DataFrame(
        [
            {
                "name": f["name"],
                "NumberOfClasses": f["NumberOfClasses"],
                "NumberOfInstances": f["NumberOfInstances"],
                "class_distribution_vector": f["class_distribution_vector"],
                "dataset": f["dataset"],
            }
            for f in filas
        ]
    )

    logger.info(
        "Cargados %d datasets de %s (modelo=%s, modalidad=%s%s).",
        len(odatasets_unique), carpeta, model, modalidad,
        f", grupo={grupo}" if grupo else "",
    )
    return odatasets_unique


if __name__ == "__main__":
    import sys
    m = sys.argv[1] if len(sys.argv) > 1 else "clip"
    mod = sys.argv[2] if len(sys.argv) > 2 else "imagen"
    df = load_all_datasets(m, mod)
    print(df[["name", "NumberOfClasses", "NumberOfInstances"]].to_string(index=False))
