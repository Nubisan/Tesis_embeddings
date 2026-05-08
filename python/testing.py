"""
testing.py
==========
Equivalente Python del archivo Testing.R.

Orquesta la ejecución de los 8 algoritmos de clustering:
  1 -> Bat
  2 -> CSCLP
  3 -> Kmedoids
  4 -> ACO
  5 -> PSO
  6 -> HCAKC
  7 -> KM-MILP
  8 -> SCK1

Modos de uso:
  python testing.py 1     # ejecuta solo el algoritmo 1
  python testing.py 9     # ejecuta los 8 en paralelo (joblib)
  python testing.py       # pregunta interactivamente

Cada algoritmo debe estar implementado en `algorithms/<modulo>.py` y exponer
una función `run(odatasets_unique)` (igual que en R, donde cada script tenía
acceso a la variable `odatasets_unique` por scoping).

Resultados (tiempo y memoria pico) se guardan en `peakRAM_log.csv`.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import os
import threading
import time
import traceback
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterator, Optional, cast

import pandas as pd
import psutil
from joblib import Parallel, delayed  # type: ignore[import-untyped]

from openml_loader import load_all_datasets

# ----------------------------------------------------------------------------
# Mapeo numalt -> módulo del algoritmo
# ----------------------------------------------------------------------------
# Cada módulo debe exponer: def run(odatasets_unique) -> Any
ALGORITHMS: dict[int, tuple[str, str]] = {
    1: ("Bat",      "algorithms.bat"),
    2: ("CSCLP",    "algorithms.csclp"),
    3: ("Kmedoids", "algorithms.kmedoids"),
    4: ("ACO",      "algorithms.aco"),
    5: ("PSO",      "algorithms.pso"),
    6: ("HCAKC",    "algorithms.hcakc"),
    7: ("KM-MILP",  "algorithms.km_milp"),
    8: ("SCK1",     "algorithms.sck1_final"),
}

PROJECT_ROOT = Path(__file__).resolve().parent
LOG_PATH     = PROJECT_ROOT / "peakRAM_log.csv"


# ----------------------------------------------------------------------------
# Medición de tiempo y RAM pico (equivalente a peakRAM de R)
# ----------------------------------------------------------------------------

@dataclass
class ResourceMeasurement:
    """Información que mide cada ejecución, equivalente a peakRAM."""
    algoritmo:          str
    elapsed_time_sec:   float
    total_ram_used_mib: float   # RAM "usada" durante la ejecución (peak - start)
    peak_ram_mib:       float   # RSS pico absoluto durante la ejecución
    fecha:              str
    estado:             str     # "OK" o "ERROR: ..."


@dataclass
class _MeasureResult:
    """Resultado de tiempo y memoria devuelto por measure_resources()."""
    elapsed:  float = 0.0
    used_mib: float = 0.0
    peak_mib: float = 0.0


class _PeakMemorySampler:
    """
    Hilo daemon que muestrea el RSS del proceso actual cada `interval`
    segundos y guarda el pico observado. Equivalente conceptual a peakRAM.
    """

    def __init__(self, interval: float = 0.05) -> None:
        self.interval = interval
        self.process  = psutil.Process(os.getpid())
        self._stop_event                       = threading.Event()
        self._thread:  Optional[threading.Thread] = None
        self.start_rss_bytes: int = 0
        self.peak_rss_bytes:  int = 0

    def _sample_loop(self) -> None:
        while not self._stop_event.is_set():
            try:
                rss = self.process.memory_info().rss
                if rss > self.peak_rss_bytes:
                    self.peak_rss_bytes = rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
            self._stop_event.wait(self.interval)

    def __enter__(self) -> "_PeakMemorySampler":
        self.start_rss_bytes = self.process.memory_info().rss
        self.peak_rss_bytes  = self.start_rss_bytes
        self._stop_event.clear()
        self._thread = threading.Thread(target=self._sample_loop, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *exc: Any) -> None:
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=1.0)

    @property
    def used_mib(self) -> float:
        """RAM "usada" durante la ejecución, en MiB (peak - start)."""
        return max(0.0, (self.peak_rss_bytes - self.start_rss_bytes) / (1024 ** 2))

    @property
    def peak_mib(self) -> float:
        """RSS pico absoluto, en MiB."""
        return self.peak_rss_bytes / (1024 ** 2)


@contextmanager
def measure_resources() -> Iterator[_MeasureResult]:
    """
    Context manager que mide tiempo y RAM. Devuelve un objeto con atributos
    `.elapsed`, `.used_mib`, `.peak_mib` una vez salido del bloque.
    """
    result  = _MeasureResult()
    sampler = _PeakMemorySampler()
    start   = time.perf_counter()
    with sampler:
        yield result
    result.elapsed  = time.perf_counter() - start
    result.used_mib = sampler.used_mib
    result.peak_mib = sampler.peak_mib


# ----------------------------------------------------------------------------
# Logging a CSV
# ----------------------------------------------------------------------------

def append_to_log(measurement: ResourceMeasurement) -> None:
    """Añade una fila a peakRAM_log.csv, creando el header si no existe."""
    write_header = not LOG_PATH.exists()
    with open(LOG_PATH, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(measurement).keys()))
        if write_header:
            writer.writeheader()
        writer.writerow(asdict(measurement))


# ----------------------------------------------------------------------------
# Ejecución de un algoritmo
# ----------------------------------------------------------------------------

def run_algoritmo(algo_id: int, odatasets_unique: pd.DataFrame) -> str:
    """
    Equivalente a `run_algoritmo_future` de R.
    Ejecuta el algoritmo `algo_id`, mide recursos y registra en el CSV.
    """
    if algo_id not in ALGORITHMS:
        return f"Opción inválida: {algo_id}"

    algo_name, module_path = ALGORITHMS[algo_id]
    print(f">>> Iniciando {algo_name} en PID {os.getpid()}", flush=True)

    # Variables siempre definidas, así no quedan unbound si falla antes del with
    estado:   str   = "OK"
    elapsed:  float = 0.0
    used_mib: float = 0.0
    peak_mib: float = 0.0

    try:
        # Importación dinámica: equivalente a source(archivo) de R
        module = importlib.import_module(module_path)

        if not hasattr(module, "run"):
            raise AttributeError(
                f"El módulo '{module_path}' debe definir "
                "una función run(odatasets_unique)"
            )

        with measure_resources() as res:
            module.run(odatasets_unique)

        # Sólo llegamos aquí si todo el bloque with se ejecutó sin excepción
        elapsed  = res.elapsed
        used_mib = res.used_mib
        peak_mib = res.peak_mib

    except Exception as e:  # pylint: disable=broad-exception-caught
        estado = f"ERROR: {e}"
        traceback.print_exc()

    # Registrar siempre, haya error o no (igual que el peakRAM de R)
    measurement = ResourceMeasurement(
        algoritmo          = algo_name,
        elapsed_time_sec   = round(elapsed, 4),
        total_ram_used_mib = round(used_mib, 2),
        peak_ram_mib       = round(peak_mib, 2),
        fecha              = datetime.now().isoformat(timespec="seconds"),
        estado             = estado,
    )
    append_to_log(measurement)

    if estado == "OK":
        return (
            f"EXITO: {algo_name} | "
            f"tiempo={elapsed:.2f}s | "
            f"used={used_mib:.1f} MiB | "
            f"peak={peak_mib:.1f} MiB"
        )
    return f"ERROR en {algo_name}: {estado}"


# ----------------------------------------------------------------------------
# Orquestador principal
# ----------------------------------------------------------------------------

def execute_test(algo_id: int) -> None:
    """
    Equivalente a `Execute_Test` de R.

    algo_id = 1..8  -> ejecuta solo ese algoritmo
    algo_id = 9     -> ejecuta los 8 en paralelo con joblib
    """
    print("Cargando datasets...", flush=True)
    odatasets_unique = load_all_datasets()

    if algo_id == 9:
        n_cores = max(1, (os.cpu_count() or 2) - 1)
        print(
            f"\n--- Ejecución paralela con joblib | {n_cores} workers ---",
            flush=True,
        )

        # backend='loky' -> procesos separados (CPU-bound).
        # joblib serializa odatasets_unique a cada worker. Para arrays grandes
        # de NumPy hace memmap automático, lo cual será útil con embeddings.
        # Usamos cast() porque joblib carece de stubs de tipos: Pylance ve un
        # tipo unión genérico (Generator | list | ...), pero sabemos que con
        # un generador como argumento devuelve list[str].
        resultados = cast(
            "list[str]",
            Parallel(n_jobs=n_cores, backend="loky", verbose=10)(
                delayed(run_algoritmo)(i, odatasets_unique)
                for i in sorted(ALGORITHMS.keys())
            ),
        )

        print("\n--- Resumen de ejecución ---")
        for r in resultados:
            print(r)
    else:
        if algo_id not in ALGORITHMS:
            raise ValueError(
                "Opción inválida. Usa 1-8 para un algoritmo o 9 para todos."
            )
        print(f"Ejecutando algoritmo {algo_id}", flush=True)
        print(run_algoritmo(algo_id, odatasets_unique))


# ----------------------------------------------------------------------------
# Entrada CLI
# ----------------------------------------------------------------------------

def _parse_args() -> int:
    parser = argparse.ArgumentParser(
        description="Ejecuta uno o todos los algoritmos de clustering."
    )
    parser.add_argument(
        "numalt",
        nargs="?",
        type=int,
        help="Número de algoritmo (1-8) o 9 para ejecutar los 8 en paralelo.",
    )
    args = parser.parse_args()

    if args.numalt is None:
        try:
            args.numalt = int(input("Introduce numalt (1-8 o 9 para todos): "))
        except (EOFError, ValueError):
            parser.error("Debes proveer un entero válido.")

    return int(args.numalt)


def _main() -> None:
    # En Windows, joblib con 'loky' requiere este guard para evitar
    # ejecuciones recursivas al arrancar workers.
    selected = _parse_args()
    execute_test(selected)


if __name__ == "__main__":
    _main()