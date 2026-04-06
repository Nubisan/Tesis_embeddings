#!/usr/bin/env python3
"""
Descarga datasets desde OpenML y crea una muestra pequeña para experimentación.

Ejemplos de uso:
    python openml_image_sampler.py --group estandar
    python openml_image_sampler.py --group baratos --max-rows 200
    python openml_image_sampler.py --datasets CIFAR-10 MNIST Fashion-MNIST

Formato esperado del JSON:
{
  "defaults": {
    "max_rows": 500,
    "output_dir": "data/openml_samples",
    "random_state": 42
  },
  "groups": {
    "estandar": ["CIFAR-10", "STL-10", "CIFAR-100"],
    "baratos": ["COIL-100", "USPS", "MNIST", "Fashion-MNIST"],
    "costosos": ["ImageNet-10", "CUB-200"]
  }
}
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import openml
import pandas as pd


# Alias útiles para nombres comunes
ALIASES = {
    "cifar10": "CIFAR_10",
    "cifar-10": "CIFAR_10",
    "cifar_10": "CIFAR_10",
    "cifar100": "CIFAR_100",
    "cifar-100": "CIFAR_100",
    "cifar_100": "CIFAR_100",
    "mnist": "mnist_784",
    "mnist_784": "mnist_784",
    "fashionmnist": "Fashion-MNIST",
    "fashion-mnist": "Fashion-MNIST",
    "fashion_mnist": "Fashion-MNIST",
    "fmnist": "Fashion-MNIST",
    "moda-mnist": "Fashion-MNIST",
    "moda_mnist": "Fashion-MNIST",
    "usps": "USPS",
    "stl10": "STL_10",
    "stl-10": "STL_10",
    "stl_10": "STL_10",
    "coil100": "coil2000",
    "coil-100": "coil2000",
    "coil_100": "coil2000",
}


def slugify(text: str) -> str:
    text = text.strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def normalize_key(text: str) -> str:
    return re.sub(r"[\s_\-]+", "", str(text).strip().lower())


def load_config(config_path: Path) -> dict[str, Any]:
    if not config_path.exists():
        raise FileNotFoundError(
            f"No se encontró el archivo de configuración: {config_path}"
        )

    with config_path.open("r", encoding="utf-8") as f:
        config = json.load(f)

    if "groups" not in config:
        config["groups"] = {}

    if "defaults" not in config:
        config["defaults"] = {}

    return config


def get_id_column(df: pd.DataFrame) -> str:
    for col in ["did", "data_id", "dataset_id"]:
        if col in df.columns:
            return col

    if df.index.name in ["did", "data_id", "dataset_id"]:
        return df.index.name

    raise ValueError(
        "No pude identificar la columna ID del dataset en la respuesta de OpenML."
    )


def get_name_column(df: pd.DataFrame) -> str:
    for col in ["name", "data_name"]:
        if col in df.columns:
            return col
    raise ValueError(
        "No pude identificar la columna nombre del dataset en la respuesta de OpenML."
    )


def resolve_dataset_id(query: str | int) -> tuple[int, str]:
    """
    Resuelve un dataset por ID o por nombre.
    Prioriza datasets activos y, si hay varias versiones, usa la más alta.
    """
    if isinstance(query, int) or (isinstance(query, str) and str(query).isdigit()):
        did = int(query)
        return did, str(did)

    raw_query = str(query).strip()
    key = normalize_key(raw_query)

    candidate_names = []
    if key in ALIASES:
        candidate_names.append(ALIASES[key])

    candidate_names.extend(
        [
            raw_query,
            raw_query.replace("-", "_"),
            raw_query.replace("_", "-"),
            raw_query.replace(" ", "_"),
            raw_query.replace(" ", "-"),
        ]
    )

    seen = set()
    candidate_names = [c for c in candidate_names if not (c in seen or seen.add(c))]

    for candidate in candidate_names:
        matches = openml.datasets.list_datasets(
            output_format="dataframe",
            status="all",
            data_name=candidate,
        )

        if matches is None or not isinstance(matches, pd.DataFrame) or matches.empty:
            continue

        matches = matches.copy()

        if "did" not in matches.columns and matches.index.name == "did":
            matches = matches.reset_index()
        elif "did" not in matches.columns:
            matches = matches.reset_index()

        id_col = get_id_column(matches)
        name_col = get_name_column(matches)

        # Intentar quedarnos con coincidencias exactas normalizadas
        norm_candidate = normalize_key(candidate)
        matches["_norm_name"] = matches[name_col].astype(str).map(normalize_key)
        exact_matches = matches[matches["_norm_name"] == norm_candidate]
        if not exact_matches.empty:
            matches = exact_matches

        # Preferir activos
        if "status" in matches.columns:
            matches["_active_rank"] = (matches["status"].astype(str) == "active").astype(int)
        else:
            matches["_active_rank"] = 0

        sort_cols = ["_active_rank"]
        ascending = [False]

        if "version" in matches.columns:
            sort_cols.append("version")
            ascending.append(False)

        matches = matches.sort_values(sort_cols, ascending=ascending)
        chosen = matches.iloc[0]

        dataset_id = int(chosen[id_col])
        dataset_name = str(chosen[name_col])
        return dataset_id, dataset_name

    raise ValueError(
        f"No pude resolver '{query}' en OpenML. "
        "Prueba con el nombre exacto en OpenML o usa directamente el data_id."
    )


def infer_target_name(dataset: openml.datasets.OpenMLDataset) -> str | None:
    target = dataset.default_target_attribute
    if target is None:
        return None

    if isinstance(target, str) and "," in target:
        return target.split(",")[0].strip()

    return str(target)


def ensure_dataframe(X: Any, attribute_names: list[str] | None = None) -> pd.DataFrame:
    if isinstance(X, pd.DataFrame):
        return X.copy()

    if attribute_names is None:
        return pd.DataFrame(X)

    return pd.DataFrame(X, columns=attribute_names)


def ensure_series(y: Any, name: str = "target") -> pd.Series | None:
    if y is None:
        return None

    if isinstance(y, pd.Series):
        out = y.copy()
        if out.name is None:
            out.name = name
        return out

    if isinstance(y, pd.DataFrame):
        if y.shape[1] == 1:
            out = y.iloc[:, 0].copy()
            if out.name is None:
                out.name = name
            return out
        raise ValueError("Se recibió un target con múltiples columnas; se esperaba una sola.")

    return pd.Series(y, name=name)


def stratified_sample_indices(
    y: pd.Series,
    n_samples: int,
    random_state: int,
) -> np.ndarray:
    """
    Muestreo estratificado robusto.
    Si n_samples < número de clases, cae a muestreo aleatorio simple.
    """
    y = y.reset_index(drop=True)
    rng = np.random.default_rng(random_state)

    groups = {}
    for class_value, indices in y.groupby(y, dropna=False).groups.items():
        idx = np.array(list(indices), dtype=int)
        rng.shuffle(idx)
        groups[class_value] = idx

    classes = list(groups.keys())
    n_classes = len(classes)

    if n_samples >= len(y):
        return np.arange(len(y), dtype=int)

    if n_samples < n_classes:
        return np.sort(rng.choice(len(y), size=n_samples, replace=False))

    selected = []
    leftovers = {}

    # Asegurar al menos una por clase
    for c in classes:
        idx = groups[c]
        selected.append(idx[0])
        leftovers[c] = idx[1:]

    remaining_budget = n_samples - n_classes
    if remaining_budget == 0:
        return np.sort(np.array(selected, dtype=int))

    leftover_counts = np.array([len(leftovers[c]) for c in classes], dtype=int)
    total_leftover = int(leftover_counts.sum())

    if total_leftover == 0:
        return np.sort(np.array(selected, dtype=int))

    desired = (leftover_counts / total_leftover) * remaining_budget
    extra = np.floor(desired).astype(int)
    extra = np.minimum(extra, leftover_counts)

    remainder = remaining_budget - int(extra.sum())
    fractional = desired - np.floor(desired)

    while remainder > 0:
        candidates = np.where(extra < leftover_counts)[0]
        if len(candidates) == 0:
            break

        best_idx = candidates[np.argmax(fractional[candidates])]
        extra[best_idx] += 1
        fractional[best_idx] = -1.0
        remainder -= 1

    for i, c in enumerate(classes):
        k = int(extra[i])
        if k > 0:
            selected.extend(leftovers[c][:k].tolist())

    return np.sort(np.array(selected, dtype=int))


def random_sample_indices(n_rows: int, n_samples: int, random_state: int) -> np.ndarray:
    rng = np.random.default_rng(random_state)
    return np.sort(rng.choice(n_rows, size=n_samples, replace=False))


def download_and_sample(
    dataset_query: str | int,
    max_rows: int,
    output_dir: Path,
    random_state: int,
    use_stratified: bool = True,
) -> dict[str, Any]:
    dataset_id, resolved_name = resolve_dataset_id(dataset_query)
    dataset = openml.datasets.get_dataset(dataset_id, download_data=True)

    target_name = infer_target_name(dataset)

    X, y, categorical_indicator, attribute_names = dataset.get_data(
        target=target_name,
        dataset_format="dataframe",
    )

    X = ensure_dataframe(X, attribute_names=attribute_names)
    y = ensure_series(y, name=target_name or "target")

    original_rows = len(X)
    sampled = False

    if max_rows is not None and max_rows > 0 and original_rows > max_rows:
        if use_stratified and y is not None and y.nunique(dropna=False) > 1:
            indices = stratified_sample_indices(
                y=y,
                n_samples=max_rows,
                random_state=random_state,
            )
        else:
            indices = random_sample_indices(
                n_rows=original_rows,
                n_samples=max_rows,
                random_state=random_state,
            )

        X = X.iloc[indices].reset_index(drop=True)
        if y is not None:
            y = y.iloc[indices].reset_index(drop=True)
        sampled = True

    dataset_name = dataset.name if getattr(dataset, "name", None) else resolved_name
    safe_name = slugify(dataset_name)
    dataset_dir = output_dir / safe_name
    dataset_dir.mkdir(parents=True, exist_ok=True)

    sample_df = X.copy()
    if y is not None:
        target_col = y.name if y.name is not None else "target"
        sample_df[target_col] = y.values

    csv_path = dataset_dir / f"{safe_name}_sample.csv"
    sample_df.to_csv(csv_path, index=False)

    parquet_path = dataset_dir / f"{safe_name}_sample.parquet"
    parquet_written = False
    try:
        sample_df.to_parquet(parquet_path, index=False)
        parquet_written = True
    except Exception:
        parquet_written = False

    metadata = {
        "query": str(dataset_query),
        "resolved_dataset_id": dataset_id,
        "resolved_dataset_name": dataset_name,
        "target_name": target_name,
        "original_rows": int(original_rows),
        "sample_rows": int(len(sample_df)),
        "n_features_without_target": int(X.shape[1]),
        "n_classes": int(y.nunique(dropna=False)) if y is not None else None,
        "sampled": sampled,
        "random_state": int(random_state),
        "categorical_feature_count": int(sum(categorical_indicator))
        if categorical_indicator is not None
        else None,
        "csv_path": str(csv_path),
        "parquet_path": str(parquet_path) if parquet_written else None,
    }

    metadata_path = dataset_dir / f"{safe_name}_metadata.json"
    with metadata_path.open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)

    return metadata


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Descarga datasets desde OpenML y crea muestras pequeñas para experimentación."
    )

    parser.add_argument(
        "--datasets",
        nargs="+",
        help="Lista de datasets por nombre o ID. Ej: CIFAR-10 MNIST 40996",
    )

    parser.add_argument(
        "--group",
        type=str,
        help="Nombre del grupo definido en el JSON. Ej: estandar",
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=Path("datasets_config.json"),
        help="Ruta al archivo de configuración JSON.",
    )

    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Número máximo de instancias a conservar por dataset. Si no se da, usa el JSON.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directorio de salida. Si no se da, usa el JSON.",
    )

    parser.add_argument(
        "--random-state",
        type=int,
        default=None,
        help="Semilla para muestreo reproducible. Si no se da, usa el JSON.",
    )

    parser.add_argument(
        "--no-stratify",
        action="store_true",
        help="Desactiva muestreo estratificado y usa muestreo aleatorio simple.",
    )

    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    config = load_config(args.config)

    defaults = config.get("defaults", {})
    groups = config.get("groups", {})

    max_rows = args.max_rows if args.max_rows is not None else defaults.get("max_rows", 500)
    output_dir = args.output_dir if args.output_dir is not None else Path(defaults.get("output_dir", "data/openml_samples"))
    random_state = args.random_state if args.random_state is not None else defaults.get("random_state", 42)

    if args.group:
        if args.group not in groups:
            available = ", ".join(groups.keys()) if groups else "(sin grupos definidos)"
            raise ValueError(
                f"El grupo '{args.group}' no existe en {args.config}. "
                f"Grupos disponibles: {available}"
            )
        dataset_list = groups[args.group]
    elif args.datasets:
        dataset_list = args.datasets
    else:
        raise ValueError("Debes usar --group o --datasets.")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []

    print(f"Configuración usada:")
    print(f"  max_rows     = {max_rows}")
    print(f"  output_dir   = {output_dir}")
    print(f"  random_state = {random_state}")
    print(f"  stratified   = {not args.no_stratify}")
    print(f"  datasets     = {dataset_list}")
    print()

    for dataset_query in dataset_list:
        try:
            metadata = download_and_sample(
                dataset_query=dataset_query,
                max_rows=max_rows,
                output_dir=output_dir,
                random_state=random_state,
                use_stratified=not args.no_stratify,
            )
            results.append(metadata)

            print(
                f"[OK] {metadata['resolved_dataset_name']} "
                f"(id={metadata['resolved_dataset_id']}) -> "
                f"{metadata['sample_rows']}/{metadata['original_rows']} filas "
                f"| features={metadata['n_features_without_target']} "
                f"| clases={metadata['n_classes']}"
            )

        except Exception as e:
            print(f"[ERROR] {dataset_query}: {e}")

    summary_name = f"summary_{args.group}.json" if args.group else "summary_manual.json"
    summary_path = output_dir / summary_name

    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print()
    print(f"Resumen guardado en: {summary_path}")


if __name__ == "__main__":
    main()