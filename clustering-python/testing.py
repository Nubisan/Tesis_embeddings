"""
testing.py
==========
Orquestador principal.

Modos de uso:

  # Un algoritmo individual con embeddings CLIP de IMAGEN (default)
  python testing.py 1 --model clip                 # 1 = BAT, imagen
  python testing.py 2 --model clip                 # 2 = CSCLP, imagen
  ...
  python testing.py 8 --model clip                 # 8 = SCK1, imagen

  # Eligiendo la MODALIDAD (segundo argumento posicional: imagen | texto)
  python testing.py 9 --model clip texto           # TODOS, CLIP, texto
  python testing.py 9 --model clip imagen          # TODOS, CLIP, imagen
  python testing.py 1 --model blip texto           # BAT, BLIP, texto

  # Filtrando solo un grupo del catalogo (A o B)
  python testing.py 9 --model qwen_vl texto --grupo A

  # Con otros modelos:
  python testing.py 9 --model blip texto
  python testing.py 9 --model qwen_vl imagen
  python testing.py 9 --model internvl texto

Si se omite la modalidad -> 'imagen'. Si se omite --model -> 'clip'.

Los CSVs de resultados se guardan en `predictions/<modelo>/<modalidad>/`,
porque testing.py setea CLUSTERING_MODEL="<modelo>" y CLUSTERING_MODALIDAD=
"<modalidad>" antes de
invocar a cada algoritmo, y `_common.get_predictions_dir()` la lee. Asi los
resultados de imagen y texto del mismo modelo NO se pisan.
"""

from __future__ import annotations

import argparse
import importlib
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Optional

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
VALID_MODALIDADES = {"imagen", "texto"}

PEAK_RAM_CSV = PROJECT_ROOT / "peakRAM_log.csv"


# ============================================================================
# Helpers
# ============================================================================

def _set_active_model(model: str, modalidad: str) -> None:
    """
    Setea las variables de entorno que los algoritmos leeran para saber donde
    guardar sus CSVs: predictions/<model>/<modalidad>/ (ver _common.get_predictions_dir).
    """
    os.environ["CLUSTERING_MODEL"] = model.lower()
    os.environ["CLUSTERING_MODALIDAD"] = modalidad.lower()


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
    algorithm_name: str, model: str, modalidad: str, elapsed: float, peak_mb: float
) -> None:
    """Append una fila al log de RAM."""
    row = {
        "timestamp":   pd.Timestamp.now().isoformat(),
        "model":       model,
        "modalidad":   modalidad,
        "algorithm":   algorithm_name,
        "elapsed_s":   round(elapsed, 4),
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
    modalidad: str,
) -> None:
    """Ejecuta un unico algoritmo."""
    if numalt not in ALGORITHMS:
        raise ValueError(f"numalt={numalt} no valido. Opciones: 1..8 o 9.")

    name, module_path = ALGORITHMS[numalt]
    print(f"\n>>> Ejecutando {name} con modelo '{model}' (modalidad '{modalidad}') ...")

    module = importlib.import_module(module_path)
    result, elapsed, peak_mb = _measure_peak_ram(module.run, odatasets_unique)
    _ = result

    _append_peak_ram_log(name, model, modalidad, elapsed, peak_mb)
    print(f">>> {name} terminado en {elapsed:.2f}s (RAM pico ~{peak_mb:.1f} MB)")


def _run_one_algorithm_in_subprocess(
    numalt: int, model: str, modalidad: str, grupo: Optional[str]
) -> tuple[int, str, float, float, str]:
    """
    Worker para joblib: ejecuta UN algoritmo de forma independiente.
    Carga sus propios datos (no serializa el DataFrame entre procesos).

    IMPORTANTE: cada worker es un proceso separado, asi que re-seteamos la
    variable de entorno aqui (no se hereda del padre con joblib/loky).
    """
    name, module_path = ALGORITHMS[numalt]
    try:
        _set_active_model(model, modalidad)  # critico: re-setear en el subproceso
        odatasets_unique = load_all_datasets(model=model, modalidad=modalidad, grupo=grupo)
        module = importlib.import_module(module_path)
        result, elapsed, peak_mb = _measure_peak_ram(
            module.run, odatasets_unique
        )
        _ = result
        return (numalt, name, elapsed, peak_mb, "OK")
    except Exception as e:  # pylint: disable=broad-exception-caught
        return (numalt, name, 0.0, 0.0, f"ERROR: {e}")


def _run_all_parallel(model: str, modalidad: str, grupo: Optional[str]) -> None:
    """Ejecuta los 8 algoritmos en paralelo con joblib."""
    print(f"\n>>> Ejecutando los 8 algoritmos en paralelo "
          f"(modelo '{model}', modalidad '{modalidad}') ...")

    n_jobs = min(8, os.cpu_count() or 4)
    print(f">>> Lanzando {n_jobs} workers (joblib + loky)...")

    results = Parallel(n_jobs=n_jobs, backend="loky", verbose=10)(
        delayed(_run_one_algorithm_in_subprocess)(numalt, model, modalidad, grupo)
        for numalt in range(1, 9)
    )

    print("\n=== Resumen ejecucion paralela ===")
    for numalt, name, elapsed, peak_mb, status in results:
        if status == "OK":
            print(
                f"  [{numalt}] {name:10s} {elapsed:7.2f}s  RAM~{peak_mb:6.1f} MB"
            )
            _append_peak_ram_log(name, model, modalidad, elapsed, peak_mb)
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
        "modalidad", type=str, nargs="?", default="imagen",
        choices=sorted(VALID_MODALIDADES),
        help="imagen | texto  (segundo argumento posicional; default: imagen)",
    )
    parser.add_argument(
        "--model", type=str, default="clip", choices=sorted(VALID_MODELS),
        help="Modelo de embeddings a usar (default: clip)",
    )
    parser.add_argument(
        "--grupo", type=str, default=None, choices=["A", "B"],
        help="Filtra solo el grupo A o B del catalogo (default: ambos)",
    )
    args = parser.parse_args()

    # CRITICO: setear el modelo+modalidad activos ANTES de invocar a cualquier
    # algoritmo para que get_predictions_dir() devuelva la ruta correcta.
    _set_active_model(args.model, args.modalidad)
    destino = f"{args.model}/{args.modalidad}"
    print(f">>> Activo: modelo '{args.model}', modalidad '{args.modalidad}'"
          f"{f', grupo {args.grupo}' if args.grupo else ''} "
          f"(resultados -> predictions/{destino}/)")

    if args.numalt == 9:
        _run_all_parallel(args.model, args.modalidad, args.grupo)
        return

    print(f"\n>>> Cargando embeddings '{args.model}' / '{args.modalidad}' ...")
    odatasets_unique = load_all_datasets(
        model=args.model, modalidad=args.modalidad, grupo=args.grupo
    )
    print(f">>> {len(odatasets_unique)} datasets cargados.")

    _run_one_algorithm(args.numalt, odatasets_unique, args.model, args.modalidad)


if __name__ == "__main__":
    main()
