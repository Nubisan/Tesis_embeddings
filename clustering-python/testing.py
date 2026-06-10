"""
testing.py
==========
Orquestador principal.

Modos de uso:

  # Ejecutar un algoritmo individual con embeddings CLIP
  python testing.py 1 --model clip          # 1 = BAT
  python testing.py 2 --model clip          # 2 = CSCLP
  ...
  python testing.py 8 --model clip          # 8 = SCK1

  # Ejecutar TODOS los algoritmos en paralelo (recomendado al final)
  python testing.py 9 --model clip          # 9 = parallel ALL

  # Con otros modelos:
  python testing.py 9 --model blip
  python testing.py 9 --model qwen_vl
  python testing.py 9 --model internvl

Si el flag --model se omite, por defecto usa clip.

Los CSVs de resultados se guardan automáticamente en `predictions/<modelo>/`,
gracias a que testing.py setea la variable de entorno CLUSTERING_MODEL antes
de invocar a cada algoritmo, y `_common.get_predictions_dir()` la lee.
"""

from __future__ import annotations

import argparse
import importlib
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd
import psutil  # type: ignore[import-untyped]
from joblib import Parallel, delayed  # type: ignore[import-untyped]

# joblib + loky preferido para paralelismo seguro con NumPy/sklearn
os.environ.setdefault("JOBLIB_TEMP_FOLDER", str(Path.cwd() / ".joblib_tmp"))

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))

# Loader de embeddings (sustituye a openml_loader)
from embeddings_loader import load_all_datasets  # noqa: E402

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(message)s")


# ============================================================================
# Mapping de algoritmos
# ============================================================================

ALGORITHMS: dict[int, tuple[str, str]] = {
    1: ("BAT",       "algorithms.bat"),
    2: ("CSCLP",     "algorithms.csclp"),
    3: ("Kmedoids",  "algorithms.kmedoids"),
    4: ("ACO",       "algorithms.aco"),
    5: ("PSO",       "algorithms.pso"),
    6: ("HCAKC",     "algorithms.hcakc"),
    7: ("KM-MILP",   "algorithms.km_milp"),
    8: ("SCK1",      "algorithms.sck1_final"),
}

VALID_MODELS = {"clip", "blip", "qwen_vl", "internvl"}

PEAK_RAM_CSV = PROJECT_ROOT / "peakRAM_log.csv"


# ============================================================================
# Helpers
# ============================================================================

def _set_active_model(model: str) -> None:
    """
    Setea la variable de entorno que los algoritmos leerán para saber
    dónde guardar sus CSVs.
    """
    os.environ["CLUSTERING_MODEL"] = model.lower()


def _measure_peak_ram(func, *args, **kwargs) -> tuple[Any, float, float]:
    """
    Ejecuta `func` y devuelve (resultado, elapsed_seconds, peak_ram_MB).
    Toma una muestra puntual de RAM (no es un sampler continuo, pero es
    suficiente para tener un orden de magnitud).
    """
    process = psutil.Process(os.getpid())
    mem_before = process.memory_info().rss / (1024 ** 2)
    start = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    mem_after = process.memory_info().rss / (1024 ** 2)
    peak_mb = max(mem_before, mem_after)
    return result, elapsed, peak_mb


def _append_peak_ram_log(
    algorithm_name: str, model: str, elapsed: float, peak_mb: float
) -> None:
    """Append una fila al log de RAM."""
    row = {
        "timestamp":  pd.Timestamp.now().isoformat(),
        "model":      model,
        "algorithm":  algorithm_name,
        "elapsed_s":  round(elapsed, 4),
        "peak_ram_MB": round(peak_mb, 2),
    }
    df_row = pd.DataFrame([row])
    if PEAK_RAM_CSV.exists():
        df_row.to_csv(PEAK_RAM_CSV, mode="a", header=False, index=False)
    else:
        df_row.to_csv(PEAK_RAM_CSV, mode="w", header=True, index=False)


# ============================================================================
# Ejecutores
# ============================================================================

def _run_one_algorithm(
    numalt: int,
    odatasets_unique: pd.DataFrame,
    model: str,
) -> None:
    """Ejecuta un único algoritmo."""
    if numalt not in ALGORITHMS:
        raise ValueError(f"numalt={numalt} no válido. Opciones: 1..8 o 9.")

    name, module_path = ALGORITHMS[numalt]
    print(f"\n>>> Ejecutando {name} con modelo '{model}' ...")

    module = importlib.import_module(module_path)
    result, elapsed, peak_mb = _measure_peak_ram(module.run, odatasets_unique)
    _ = result

    _append_peak_ram_log(name, model, elapsed, peak_mb)
    print(f">>> {name} terminado en {elapsed:.2f}s (RAM pico ~{peak_mb:.1f} MB)")


def _run_one_algorithm_in_subprocess(
    numalt: int, model: str
) -> tuple[int, str, float, float, str]:
    """
    Worker para joblib: ejecuta UN algoritmo de forma independiente.
    Carga sus propios datos (no serializa el DataFrame entre procesos —
    eso sería más lento por la matriz de embeddings).

    IMPORTANTE: cada worker es un proceso separado, así que tenemos que
    re-setear la variable de entorno aquí (no se hereda del padre porque
    joblib/loky no copia el environment dinámicamente).
    """
    name, module_path = ALGORITHMS[numalt]
    try:
        _set_active_model(model)  # crítico: re-setear en el subproceso
        odatasets_unique = load_all_datasets(model=model)
        module = importlib.import_module(module_path)
        result, elapsed, peak_mb = _measure_peak_ram(
            module.run, odatasets_unique
        )
        _ = result
        return (numalt, name, elapsed, peak_mb, "OK")
    except Exception as e:  # pylint: disable=broad-exception-caught
        return (numalt, name, 0.0, 0.0, f"ERROR: {e}")


def _run_all_parallel(model: str) -> None:
    """Ejecuta los 8 algoritmos en paralelo con joblib."""
    print(f"\n>>> Ejecutando los 8 algoritmos en paralelo con modelo '{model}' ...")

    n_jobs = min(8, os.cpu_count() or 4)
    print(f">>> Lanzando {n_jobs} workers (joblib + loky)...")

    results = Parallel(n_jobs=n_jobs, backend="loky", verbose=10)(
        delayed(_run_one_algorithm_in_subprocess)(numalt, model)
        for numalt in range(1, 9)
    )

    print("\n=== Resumen ejecución paralela ===")
    for numalt, name, elapsed, peak_mb, status in results:
        if status == "OK":
            print(
                f"  [{numalt}] {name:10s} {elapsed:7.2f}s  RAM~{peak_mb:6.1f} MB"
            )
            _append_peak_ram_log(name, model, elapsed, peak_mb)
        else:
            print(f"  [{numalt}] {name:10s} {status}")


# ============================================================================
# Main
# ============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Orquestador de algoritmos de clustering con embeddings",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "numalt", type=int,
        help="1..8 = algoritmo individual; 9 = todos en paralelo",
    )
    parser.add_argument(
        "--model", type=str, default="clip", choices=sorted(VALID_MODELS),
        help="Modelo de embeddings a usar (default: clip)",
    )
    args = parser.parse_args()

    # CRÍTICO: setear el modelo activo ANTES de invocar a cualquier algoritmo
    # para que get_predictions_dir() en _common.py devuelva la ruta correcta.
    _set_active_model(args.model)
    print(f">>> Modelo activo: '{args.model}' (resultados → predictions/{args.model}/)")

    if args.numalt == 9:
        _run_all_parallel(args.model)
        return

    print(f"\n>>> Cargando embeddings con modelo '{args.model}' ...")
    odatasets_unique = load_all_datasets(model=args.model)
    print(f">>> {len(odatasets_unique)} datasets cargados.")

    _run_one_algorithm(args.numalt, odatasets_unique, args.model)


if __name__ == "__main__":
    main()
