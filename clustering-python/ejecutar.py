"""
run_pendientes.py
=================
Orquestador para terminar las corridas pendientes (grupo B) repartidas entre
varias maquinas, SIN tocar ninguno de los 8 algoritmos ni testing.py.

Idea clave
----------
Cada algoritmo, al terminar su run(), escribe TODO su pred_<ALGO>.csv de una
sola vez. Si una corrida larga se mata a la mitad, se pierde todo lo avanzado.
Para evitarlo, este runner corre **un dataset a la vez** aprovechando la
variable de entorno EMB_ROOT que ya lee embeddings_loader:

  1. Crea un EMB_ROOT temporal con UN solo .npz (symlink al real).
  2. Invoca  `python testing.py <numalt> <modalidad> --model <model>`  con timeout.
     - numalt 9  -> los 8 algoritmos en paralelo sobre ese unico dataset.
     - numalt k  -> solo el algoritmo k (para reintentar los pesados).
  3. Cada algoritmo deja su pred_<ALGO>.csv (1 fila) en predictions/<model>/<mod>/.
     El runner lee esa fila y la ACUMULA en .runner_state/results/... (upsert por
     nombre de dataset). Esa acumulacion es la fuente de verdad.

Asi:
  * Checkpointing real: si el proceso muere, lo ya guardado queda.
  * Timeout por dataset: un MILP/PAM colgado en cifar100 se mata y no bloquea al resto.
  * Resumible: al relanzar, salta lo ya acumulado.
  * Reparto por maquina: cada maquina corre un subconjunto (--models / --modalidades).

Fases
-----
  --phase 1  : numalt 9 sobre TODOS los datasets del bucket, timeout corto (T1).
               Termina de golpe lo barato/mediano en todos los datasets.
  --phase 2  : reintenta SOLO los (algoritmo,dataset) pesados que quedaron
               pendientes, uno por uno, con timeout largo (T2).
  --merge    : fusiona grupo A (backup) + grupo B (acumulado) en los
               predictions/<model>/<mod>/pred_<ALGO>.csv finales.

Ejemplos
--------
  # Maquina 1 (rapida): fase 1 de imagen para clip y blip
  python run_pendientes.py --phase 1 --models clip blip --modalidades imagen \
      --emb-root ./embeddings --t1 1200

  # Maquina 2: fase 2 (pesados) de imagen para clip y blip
  python run_pendientes.py --phase 2 --models clip blip --modalidades imagen \
      --emb-root ./embeddings --t2 7200

  # Al final, en cada maquina (o donde se junten los .runner_state):
  python run_pendientes.py --merge --models clip blip --modalidades imagen
"""
from __future__ import annotations

import argparse
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent
STATE_DIR = PROJECT_ROOT / ".runner_state"
RESULTS_DIR = STATE_DIR / "results"          # acumulador (fuente de verdad, grupo B)
TMP_ROOT = STATE_DIR / "tmp_emb"             # EMB_ROOT temporales de 1 dataset
LOG_CSV = STATE_DIR / "runner_log.csv"
BACKUP_DIR = PROJECT_ROOT / "predictions_backup_grupoA"
PREDICTIONS_DIR = PROJECT_ROOT / "predictions"

# numalt -> nombre de archivo pred_<ALGO>.csv que produce cada algoritmo
ALGO_FILE = {
    1: "pred_BAT.csv",
    2: "pred_CSCLP.csv",
    3: "pred_KMEDOIDS.csv",
    4: "pred_ACO.csv",
    5: "pred_PSO.csv",
    6: "pred_HCAKC.csv",
    7: "pred_KM-MILP.csv",
    8: "pred_SCK1.csv",
}
ALL_ALGOS = list(ALGO_FILE.keys())
# Orden de reintento en fase 2: baratos primero (mop-up rapido), pesados al
# final (donde el timeout largo importa). Es red de seguridad de TODO lo que
# falte, no solo de los pesados.
PHASE2_ORDER = [8, 5, 6, 1, 7, 4, 2, 3]  # SCK1,PSO,HCAKC,BAT,KM-MILP,ACO,CSCLP,Kmedoids

MODELS = ["clip", "blip", "qwen_vl", "internvl"]
MODALIDADES = ["imagen", "texto"]


# ---------------------------------------------------------------------------
# Utilidades de estado
# ---------------------------------------------------------------------------

def _acc_path(model: str, modalidad: str, numalt: int) -> Path:
    return RESULTS_DIR / model / modalidad / ALGO_FILE[numalt]


def _done_datasets(model: str, modalidad: str, numalt: int) -> set[str]:
    """Nombres de dataset ya acumulados para (model, modalidad, algoritmo)."""
    p = _acc_path(model, modalidad, numalt)
    if not p.exists():
        return set()
    try:
        return set(pd.read_csv(p)["name"].astype(str))
    except Exception:
        return set()


def _upsert_row(model: str, modalidad: str, numalt: int, row: pd.DataFrame) -> None:
    """Inserta/actualiza (por 'name') las filas de 'row' en el acumulador."""
    p = _acc_path(model, modalidad, numalt)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        base = pd.read_csv(p)
        base = base[~base["name"].astype(str).isin(row["name"].astype(str))]
        out = pd.concat([base, row], ignore_index=True)
    else:
        out = row
    out.to_csv(p, index=False)


def _log(**kw) -> None:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    kw = {"timestamp": pd.Timestamp.now().isoformat(), **kw}
    df = pd.DataFrame([kw])
    df.to_csv(LOG_CSV, mode="a", header=not LOG_CSV.exists(), index=False)


def _list_datasets(emb_root: Path, model: str, modalidad: str,
                   grupo: str | None = None) -> list[Path]:
    folder = emb_root / model / modalidad
    if not folder.is_dir():
        return []
    files = sorted(folder.glob("*.npz"))
    if grupo is not None:
        files = [f for f in files if _dataset_grupo(f) == str(grupo)]
    return files


def _dataset_size(npz_path: Path) -> int:
    """N de instancias (para ordenar chico->grande). Robusto a fallos."""
    try:
        import numpy as np
        with np.load(npz_path, allow_pickle=True) as z:
            return int(z["y"].shape[0])
    except Exception:
        return 0


def _dataset_grupo(npz_path: Path) -> str:
    """Lee el campo 'grupo' del .npz ('' si no existe)."""
    try:
        import numpy as np
        with np.load(npz_path, allow_pickle=True) as z:
            return str(z["grupo"]) if "grupo" in z else ""
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# Ejecucion de un dataset (con timeout que mata TODO el arbol de procesos)
# ---------------------------------------------------------------------------

def _run_single(model: str, modalidad: str, numalt: int,
                npz_path: Path, timeout_s: float) -> tuple[str, float]:
    """
    Corre testing.py con un EMB_ROOT que contiene solo `npz_path`.
    numalt=9 -> los 8 algoritmos; numalt k -> solo ese.
    Devuelve (status, elapsed). status in {"OK","TIMEOUT","ERROR"}.
    """
    name = npz_path.stem
    tmp = TMP_ROOT / f"{model}__{modalidad}__{numalt}__{name}" / model / modalidad
    tmp.mkdir(parents=True, exist_ok=True)
    link = tmp / npz_path.name
    if link.exists() or link.is_symlink():
        link.unlink()
    try:
        link.symlink_to(npz_path.resolve())
    except OSError:
        # Si el FS no soporta symlink (algunos Windows), copiar.
        import shutil
        shutil.copy2(npz_path, link)

    env = dict(os.environ)
    env["EMB_ROOT"] = str((TMP_ROOT / f"{model}__{modalidad}__{numalt}__{name}").resolve())

    cmd = [sys.executable, "testing.py", str(numalt), modalidad, "--model", model]
    is_windows = (os.name == "nt")
    # Crear un grupo de procesos propio para poder matar joblib + TODOS sus hijos.
    if is_windows:
        popen_kw = {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP}
    else:
        popen_kw = {"start_new_session": True}

    start = time.perf_counter()
    proc = subprocess.Popen(cmd, cwd=str(PROJECT_ROOT), env=env,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                            **popen_kw)
    try:
        proc.wait(timeout=timeout_s)
        status = "OK"
    except subprocess.TimeoutExpired:
        _kill_tree(proc, is_windows)
        proc.wait()
        status = "TIMEOUT"
    except Exception:
        status = "ERROR"
    elapsed = time.perf_counter() - start
    return status, elapsed


def _kill_tree(proc: subprocess.Popen, is_windows: bool) -> None:
    """Mata el proceso y TODO su arbol (joblib/loky lanza procesos hijos)."""
    try:
        if is_windows:
            # /T mata el arbol; /F fuerza. Evita workers de loky huerfanos.
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except (ProcessLookupError, OSError):
        pass


def _harvest(model: str, modalidad: str, name: str, numalts: list[int]) -> list[int]:
    """
    Lee los pred_<ALGO>.csv scratch en predictions/<model>/<mod>/ y acumula la
    fila de `name` para cada algoritmo que la haya producido.
    Devuelve la lista de numalts efectivamente cosechados.
    """
    harvested: list[int] = []
    scratch = PREDICTIONS_DIR / model / modalidad
    for k in numalts:
        f = scratch / ALGO_FILE[k]
        if not f.exists():
            continue
        try:
            df = pd.read_csv(f)
        except Exception:
            continue
        row = df[df["name"].astype(str) == str(name)]
        if len(row) >= 1:
            _upsert_row(model, modalidad, k, row)
            harvested.append(k)
    return harvested


# ---------------------------------------------------------------------------
# Fases
# ---------------------------------------------------------------------------

def phase1(models, modalidades, emb_root: Path, t1: float, grupo=None) -> None:
    """numalt 9 sobre cada dataset del bucket; timeout corto por dataset."""
    for model in models:
        for modalidad in modalidades:
            datasets = _list_datasets(emb_root, model, modalidad, grupo)
            if not datasets:
                print(f"[!] Sin .npz (grupo={grupo}) en {emb_root}/{model}/{modalidad} -- salto.")
                continue
            # chico -> grande: los riesgosos (grandes) al final
            datasets.sort(key=_dataset_size)
            for npz in datasets:
                name = npz.stem
                pend = [k for k in ALL_ALGOS
                        if name not in _done_datasets(model, modalidad, k)]
                if not pend:
                    print(f"  [skip] {model}/{modalidad}/{name} (completo)")
                    continue
                print(f"  [F1] {model}/{modalidad}/{name} (n={_dataset_size(npz)}) "
                      f"faltan {len(pend)} algos ...", flush=True)
                status, el = _run_single(model, modalidad, 9, npz, t1)
                got = _harvest(model, modalidad, name, ALL_ALGOS)
                miss = [k for k in ALL_ALGOS
                        if name not in _done_datasets(model, modalidad, k)]
                _log(phase=1, model=model, modalidad=modalidad, dataset=name,
                     status=status, elapsed_s=round(el, 1),
                     harvested=len(got), missing=len(miss),
                     missing_algos=";".join(ALGO_FILE[k].replace('pred_', '').replace('.csv', '') for k in miss))
                tag = "OK" if not miss else f"faltan {[ALGO_FILE[k][5:-4] for k in miss]}"
                print(f"       {status} {el:.0f}s -> cosechados {len(got)}; {tag}")


def phase2(models, modalidades, emb_root: Path, t2: float,
           order=PHASE2_ORDER, grupo=None) -> None:
    """Reintenta uno-por-uno CADA (algoritmo, dataset) aun pendiente.
    Baratos primero; pesados al final con el timeout largo t2."""
    for model in models:
        for modalidad in modalidades:
            datasets = _list_datasets(emb_root, model, modalidad, grupo)
            datasets.sort(key=_dataset_size)          # chico -> grande
            for npz in datasets:
                name = npz.stem
                for k in order:
                    if name in _done_datasets(model, modalidad, k):
                        continue
                    algo = ALGO_FILE[k][5:-4]
                    print(f"  [F2] {model}/{modalidad}/{name} (n={_dataset_size(npz)}) "
                          f"<- {algo} (timeout {t2:.0f}s) ...", flush=True)
                    status, el = _run_single(model, modalidad, k, npz, t2)
                    got = _harvest(model, modalidad, name, [k])
                    ok = bool(got)
                    _log(phase=2, model=model, modalidad=modalidad, dataset=name,
                         algo=algo, status=status if ok else "NO_RESULT",
                         elapsed_s=round(el, 1))
                    print(f"       {'HECHO' if ok else 'sin resultado ('+status+')'} {el:.0f}s")


def merge(models, modalidades) -> None:
    """
    Final = grupo A (backup) UNION grupo B (acumulado), dedup por 'name'
    (grupo B pisa si hay choque). Escribe predictions/<model>/<mod>/pred_<ALGO>.csv.
    """
    for model in models:
        for modalidad in modalidades:
            for k in ALL_ALGOS:
                b = _acc_path(model, modalidad, k)          # grupo B nuevo
                a = BACKUP_DIR / model / modalidad / ALGO_FILE[k]   # grupo A backup
                frames = []
                if a.exists():
                    frames.append(pd.read_csv(a))
                if b.exists():
                    frames.append(pd.read_csv(b))
                if not frames:
                    continue
                allrows = pd.concat(frames, ignore_index=True)
                # grupo B (ultimo) gana ante choque de 'name'
                allrows = allrows.drop_duplicates(subset="name", keep="last")
                out = PREDICTIONS_DIR / model / modalidad / ALGO_FILE[k]
                out.parent.mkdir(parents=True, exist_ok=True)
                allrows.to_csv(out, index=False)
                print(f"  [merge] {out}  ({len(allrows)} datasets)")


def status_report(models, modalidades, emb_root: Path, grupo=None) -> None:
    """Tabla de cuantos algos estan listos por dataset."""
    for model in models:
        for modalidad in modalidades:
            datasets = _list_datasets(emb_root, model, modalidad, grupo)
            if not datasets:
                continue
            print(f"\n== {model}/{modalidad} ==")
            for npz in sorted(datasets, key=_dataset_size):
                name = npz.stem
                done = [ALGO_FILE[k][5:-4] for k in ALL_ALGOS
                        if name in _done_datasets(model, modalidad, k)]
                print(f"  {name:16s} n={_dataset_size(npz):6d}  {len(done)}/8  "
                      f"faltan: {[ALGO_FILE[k][5:-4] for k in ALL_ALGOS if ALGO_FILE[k][5:-4] not in done]}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description="Runner de corridas pendientes (grupo B)")
    ap.add_argument("--phase", type=int, choices=[1, 2], default=None)
    ap.add_argument("--merge", action="store_true")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--models", nargs="+", default=MODELS, choices=MODELS)
    ap.add_argument("--modalidades", nargs="+", default=MODALIDADES, choices=MODALIDADES)
    ap.add_argument("--emb-root", type=str, default=str(PROJECT_ROOT / "embeddings"))
    ap.add_argument("--t1", type=float, default=1200.0, help="timeout/dataset fase 1 (s)")
    ap.add_argument("--t2", type=float, default=7200.0, help="timeout/algo-dataset fase 2 (s)")
    ap.add_argument("--grupo", type=str, default=None, choices=["A", "B"],
                    help="procesa solo datasets de este grupo (lee el campo grupo del .npz)")
    ap.add_argument("--km-timelimit", type=float, default=None,
                    help="time_limit (s) por resolución MILP de KM-MILP; usa la incumbente si se alcanza")
    ap.add_argument("--backup", action="store_true",
                    help="respalda predictions/ -> predictions_backup_grupoA/ y sale")
    args = ap.parse_args()

    emb_root = Path(args.emb_root).resolve()

    # KM-MILP: propagar el time_limit a los subprocesos via env var.
    if args.km_timelimit:
        os.environ["KM_MILP_TIME_LIMIT"] = str(args.km_timelimit)

    if args.backup:
        import shutil
        if PREDICTIONS_DIR.exists():
            if BACKUP_DIR.exists():
                print(f"[!] {BACKUP_DIR} ya existe; no lo piso.")
            else:
                shutil.copytree(PREDICTIONS_DIR, BACKUP_DIR)
                print(f"[backup] {PREDICTIONS_DIR} -> {BACKUP_DIR}")
        else:
            print("[!] No hay predictions/ que respaldar.")
        return

    if args.status:
        status_report(args.models, args.modalidades, emb_root, args.grupo)
        return
    if args.merge:
        merge(args.models, args.modalidades)
        return
    if args.phase == 1:
        phase1(args.models, args.modalidades, emb_root, args.t1, args.grupo)
        return
    if args.phase == 2:
        phase2(args.models, args.modalidades, emb_root, args.t2, grupo=args.grupo)
        return
    ap.print_help()


if __name__ == "__main__":
    main()
