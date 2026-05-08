"""
openml_loader.py
================
Equivalente Python del archivo Openml.R.

Descarga, filtra y cachea datasets de OpenML siguiendo exactamente la misma
lógica que el script de R, para asegurar resultados idénticos:

  - Filtros sobre la lista: 4-6 features, 150-250 instancias, 2-3 clases
  - distinct por (NumberOfClasses, NumberOfInstances)
  - distinct por nombre case-insensitive
  - NumberOfSymbolicFeatures no nulo y <= 1
  - NumberOfMissingValues no nulo e igual a 0
  - Validación de num. de clases real vs reportado
  - Filtro final por firma digital de distribución de clases

Se cachea:
  - La lista filtrada de OpenML en datasets_local/oml_list_cache.pkl
  - Cada dataset descargado en datasets_local/dataset_<id>.pkl
"""

from __future__ import annotations

import logging
import pickle
import time
import warnings
from pathlib import Path
from typing import Any, Optional

import pandas as pd
import openml

# ----------------------------------------------------------------------------
# Configuración global
# ----------------------------------------------------------------------------

logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

# Carpeta de caché local: se crea junto al archivo del loader, no al cwd.
# Así evitamos el bug del R original donde la ruta usaba ".." y "datasets_local"
# de forma inconsistente.
PROJECT_ROOT = Path(__file__).resolve().parent
DATASETS_DIR = PROJECT_ROOT / "datasets_local"
DATASETS_DIR.mkdir(parents=True, exist_ok=True)

LOCAL_LIST_PATH = DATASETS_DIR / "oml_list_cache.pkl"


# ----------------------------------------------------------------------------
# Filtrado de metadatos
# ----------------------------------------------------------------------------

def filter_unique_datasets_metadata(metadata_df: pd.DataFrame) -> pd.DataFrame:
    """
    Equivalente a `filter_unique_datasets_metadata` de R.

    Aplica en orden:
      1. distinct por (NumberOfClasses, NumberOfInstances)  -> primer match
      2. distinct por nombre (case-insensitive)             -> primer match
      3. Filtro NumberOfSymbolicFeatures no nulo y <= 1
      4. Filtro NumberOfMissingValues no nulo e igual a 0
    """
    out = metadata_df.copy()

    # Paso 1: distinct(NumberOfClasses, NumberOfInstances, .keep_all = TRUE)
    out = out.drop_duplicates(
        subset=["NumberOfClasses", "NumberOfInstances"], keep="first"
    )

    # Paso 2: distinct por nombre case-insensitive
    out["__name_lower"] = out["name"].str.lower()
    out = out.drop_duplicates(subset=["__name_lower"], keep="first")
    out = out.drop(columns=["__name_lower"])

    # Paso 3: NumberOfSymbolicFeatures no nulo y <= 1
    out = out[out["NumberOfSymbolicFeatures"].notna()]
    out = out[out["NumberOfSymbolicFeatures"] <= 1]

    # Paso 4: NumberOfMissingValues no nulo y == 0
    out = out[out["NumberOfMissingValues"].notna()]
    out = out[out["NumberOfMissingValues"] == 0]

    return out.reset_index(drop=True)


# ----------------------------------------------------------------------------
# Descarga de un dataset individual
# ----------------------------------------------------------------------------

def load_dataset(dataset_id: int) -> Optional[pd.DataFrame]:
    """
    Carga un dataset desde caché local; si no existe, lo descarga de OpenML.
    Equivalente a `load_dataset(id)` de R.

    Devuelve un DataFrame con todas las columnas (features + target) o None si
    hubo un error.
    """
    local_path = DATASETS_DIR / f"dataset_{dataset_id}.pkl"

    if local_path.exists():
        return pd.read_pickle(local_path)

    try:
        odata = openml.datasets.get_dataset(
            dataset_id,
            download_data=True,
            download_qualities=False,
            download_features_meta_data=False,
        )

        # get_data devuelve (X, y, categorical_indicator, attribute_names).
        # Pasando target=None obtenemos todo el dataframe junto, igual que en R.
        # Aseguramos con assert que sea DataFrame para que el type checker
        # no se confunda (puede ser ndarray o csr_matrix en otros formatos).
        data_df, _, _, _ = odata.get_data(dataset_format="dataframe", target=None)
        assert isinstance(data_df, pd.DataFrame), (
            f"Se esperaba un DataFrame para el dataset {dataset_id}, "
            f"se recibió {type(data_df).__name__}"
        )

        # Manejar columnas duplicadas tal como el R original
        if data_df.columns.duplicated().any():
            warnings.warn(
                f"Dataset ID {dataset_id}: columnas duplicadas. Renombrando..."
            )
            new_cols: list[str] = []
            counts: dict[str, int] = {}
            for col in data_df.columns:
                col_str = str(col)
                if col_str not in counts:
                    counts[col_str] = 0
                    new_cols.append(col_str)
                else:
                    counts[col_str] += 1
                    new_cols.append(f"{col_str}.{counts[col_str]}")
            data_df.columns = new_cols

        data_df.to_pickle(local_path)
        logger.info("Dataset ID %s descargado y guardado localmente.", dataset_id)
        return data_df

    except (openml.exceptions.OpenMLServerException, OSError, ValueError) as e:
        warnings.warn(f"Error cargando dataset ID {dataset_id}: {e}")
        return None


# ----------------------------------------------------------------------------
# Listado filtrado de OpenML (con caché)
# ----------------------------------------------------------------------------

def list_oml_datasets_filtered() -> pd.DataFrame:
    """
    Equivalente a la sección de `list_oml_data` de R.

    Devuelve un DataFrame con los datasets de OpenML que cumplen:
        4 <= NumberOfFeatures  <= 6
      150 <= NumberOfInstances <= 250
        2 <= NumberOfClasses   <= 3

    Usa caché local en `oml_list_cache.pkl`.
    """
    if LOCAL_LIST_PATH.exists():
        logger.info("Cargando lista de datasets desde caché local...")
        with open(LOCAL_LIST_PATH, "rb") as f:
            cached: pd.DataFrame = pickle.load(f)
        return cached

    logger.info("Consultando API de OpenML para listar datasets...")
    try:
        # En versiones recientes de openml-python list_datasets puede devolver
        # dict o DataFrame. Forzamos DataFrame con un fallback robusto.
        try:
            all_datasets = openml.datasets.list_datasets(output_format="dataframe")
        except TypeError:
            raw = openml.datasets.list_datasets()
            if isinstance(raw, dict):
                all_datasets = pd.DataFrame.from_dict(raw, orient="index")
            else:
                all_datasets = raw

        assert isinstance(all_datasets, pd.DataFrame), (
            "list_datasets no devolvió un DataFrame válido"
        )

        # Convertimos las columnas numéricas (algunas vienen como object)
        for col in [
            "NumberOfFeatures",
            "NumberOfInstances",
            "NumberOfClasses",
            "NumberOfMissingValues",
            "NumberOfSymbolicFeatures",
        ]:
            if col in all_datasets.columns:
                all_datasets[col] = pd.to_numeric(
                    all_datasets[col], errors="coerce"
                )

        # Filtros equivalentes a number_features=c(4,6), etc.
        mask = (
            (all_datasets["NumberOfFeatures"] >= 4)
            & (all_datasets["NumberOfFeatures"] <= 6)
            & (all_datasets["NumberOfInstances"] >= 150)
            & (all_datasets["NumberOfInstances"] <= 250)
            & (all_datasets["NumberOfClasses"] >= 2)
            & (all_datasets["NumberOfClasses"] <= 3)
        )
        filtered = all_datasets[mask].copy()

        # Aseguramos que exista una columna 'data_id' para alinearnos con R
        if "data_id" not in filtered.columns and "did" in filtered.columns:
            filtered["data_id"] = filtered["did"]
        elif "did" not in filtered.columns and "data_id" in filtered.columns:
            filtered["did"] = filtered["data_id"]

        with open(LOCAL_LIST_PATH, "wb") as f:
            pickle.dump(filtered, f)
        logger.info("Lista guardada en caché local.")
        return filtered

    except (openml.exceptions.OpenMLServerException, OSError) as e:
        raise RuntimeError(
            "Error crítico: no se pudo conectar con OpenML para obtener la "
            f"lista y no hay caché local. Detalle: {e}"
        ) from e


# ----------------------------------------------------------------------------
# Distribución de clases por dataset
# ----------------------------------------------------------------------------

def get_class_distributions(
    dataset_id: int, expected_classes: int
) -> dict[str, Any]:
    """
    Equivalente a `get_class_distributions` de R.

    Carga el dataset, busca la columna 'class' (o usa la última columna) como
    target, calcula la distribución de clases y valida que coincida con el
    número esperado.

    Devuelve un dict con:
      - text:   "valor1:n1; valor2:n2; ..."
      - vector: lista de enteros con los conteos
      - data:   el DataFrame completo (o None en caso de error)
    """
    dataset = load_dataset(dataset_id)

    if dataset is None:
        return {"text": "Error", "vector": None, "data": None}

    if "class" in dataset.columns:
        target = dataset["class"]
    else:
        target = dataset.iloc[:, -1]

    # value_counts ordenado por nombre del nivel, tal como `table()` de R
    class_dist = target.value_counts().sort_index()
    num_classes_found = len(class_dist)

    if num_classes_found != expected_classes:
        warnings.warn(
            f"Dataset ID {dataset_id}: discrepancia en número de clases "
            f"(esperado={expected_classes}, encontrado={num_classes_found})"
        )
        return {"text": "Error", "vector": None, "data": None}

    text = "; ".join(
        f"{name}:{int(count)}" for name, count in class_dist.items()
    )
    vector = [int(c) for c in class_dist.values]

    return {"text": text, "vector": vector, "data": dataset}


# ----------------------------------------------------------------------------
# Punto de entrada principal
# ----------------------------------------------------------------------------

def load_all_datasets() -> pd.DataFrame:
    """
    Equivalente a la ejecución completa de Openml.R.

    Devuelve `odatasets_unique`: un DataFrame con metadatos y los datasets
    cargados en una columna 'dataset' (lista de DataFrames).
    """
    odatasets = list_oml_datasets_filtered()
    odatasets_unique = filter_unique_datasets_metadata(odatasets)

    # Iteramos sobre cada (data_id, NumberOfClasses)
    start = time.time()
    class_distributions: list[dict[str, Any]] = []
    for did, n_cls in zip(
        odatasets_unique["data_id"].astype(int),
        odatasets_unique["NumberOfClasses"].astype(int),
    ):
        class_distributions.append(get_class_distributions(int(did), int(n_cls)))
    elapsed = time.time() - start
    logger.info("Total loading execution time: %.2fs", elapsed)

    # Filtramos los inválidos (data is None)
    valid_mask: list[bool] = [cd["data"] is not None for cd in class_distributions]
    odatasets_unique = odatasets_unique[valid_mask].reset_index(drop=True)
    valid: list[dict[str, Any]] = [
        cd for cd, ok in zip(class_distributions, valid_mask) if ok
    ]

    odatasets_unique["class_distribution"] = [cd["text"] for cd in valid]
    odatasets_unique["class_distribution_vector"] = [cd["vector"] for cd in valid]
    odatasets_unique["dataset"] = [cd["data"] for cd in valid]

    # Firma digital: ordenamos los conteos y los unimos por '-'
    def _signature(v: Optional[list[int]]) -> Optional[str]:
        if v is None:
            return None
        return "-".join(str(x) for x in sorted(v))

    odatasets_unique["class_dist_signature"] = (
        odatasets_unique["class_distribution_vector"].apply(_signature)
    )

    odatasets_unique = odatasets_unique.drop_duplicates(
        subset=["NumberOfClasses", "NumberOfInstances", "class_dist_signature"],
        keep="first",
    ).reset_index(drop=True)

    logger.info(
        "Datasets finales tras firma digital: %d", len(odatasets_unique)
    )
    return odatasets_unique


# ----------------------------------------------------------------------------
# Ejecución directa (para probar el loader sin pasar por testing.py)
# ----------------------------------------------------------------------------

def _main() -> None:
    result_df = load_all_datasets()
    cols_to_show = [
        c for c in [
            "data_id", "name", "NumberOfClasses",
            "NumberOfInstances", "class_distribution",
        ] if c in result_df.columns
    ]
    print("\n=== Datasets cargados ===")
    print(result_df[cols_to_show].to_string(index=False))


if __name__ == "__main__":
    _main()