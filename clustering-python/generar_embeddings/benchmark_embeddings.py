# -*- coding: utf-8 -*-
"""
benchmark_embeddings.py
=======================
Genera los embeddings en CPU midiendo, POR DATASET:
  - tiempo de carga del dataset (t_datos_s)         -> incluye descarga la 1a vez
  - tiempo puro de embedding / forward passes (t_embed_s)  <- el numero para la tesis
  - segundos por item y items por segundo (throughput, comparable entre modelos)
  - pico de RAM del proceso durante el embedding (ram_pico_MB)

Escribe un CSV incremental (una fila por dataset, se guarda al terminar cada uno)
y un JSON con la huella del entorno (CPU/RAM/versiones) para reproducibilidad.

Es REANUDABLE: si el CSV ya tiene un dataset con estado 'ok', lo salta. Puedes
cortar y volver a lanzar, o repartir datasets entre tus 3 maquinas.

EJEMPLOS
--------
# Modelos ligeros: corren completos en CPU en un rato
python benchmark_embeddings.py --modelo clip
python benchmark_embeddings.py --modelo blip

# Modelos pesados (2B/1B params): en CPU son LENTOS. Recomendado medir con tope
# de N para obtener throughput representativo sin esperar dias:
python benchmark_embeddings.py --modelo qwen_vl   --max-n 500
python benchmark_embeddings.py --modelo internvl  --max-n 500 --dtype bfloat16

# Solo una modalidad / un grupo / unos datasets:
python benchmark_embeddings.py --modelo clip --modalidad texto --grupo B
python benchmark_embeddings.py --modelo clip --solo mnist,cifar10

# Registrar el nombre exacto del CPU en la tesis:
python benchmark_embeddings.py --modelo clip --cpu-nombre "AMD Ryzen 7 5700U with Radeon Graphics"
"""

import argparse
import os
import sys
import time
import ssl

ssl._create_default_https_context = ssl._create_unverified_context
os.environ["TOKENIZERS_PARALLELISM"] = "false"
os.environ["HF_DATASETS_DISABLE_PROGRESS_BARS"] = "1"
os.environ["TQDM_DISABLE"] = "1"

try:
    import tqdm
    tqdm.tqdm.monitor_interval = 0
except Exception:
    pass


def parse_args():
    p = argparse.ArgumentParser(description='Benchmark de generacion de embeddings en CPU.')
    p.add_argument('--modelo', required=True, choices=['clip', 'blip', 'qwen_vl', 'internvl'])
    p.add_argument('--device', default='auto', choices=['auto', 'cpu', 'cuda'],
                   help='auto = cuda si hay GPU, si no cpu')
    p.add_argument('--modalidad', default='ambas', choices=['imagen', 'texto', 'ambas'])
    p.add_argument('--grupo', default='todos', choices=['A', 'B', 'todos'])
    p.add_argument('--solo', default='', help='lista separada por comas de datasets a correr')
    p.add_argument('--max-n', type=int, default=0,
                   help='tope de instancias por dataset (submuestreo estratificado, seed=123). 0=sin tope')
    p.add_argument('--dtype', default=None, choices=['float32', 'bfloat16'],
                   help='precision en CPU (por defecto float32; bfloat16 ayuda a Qwen a entrar en 16 GB)')
    p.add_argument('--hilos', type=int, default=0, help='hilos de CPU (0 = nucleos fisicos)')
    p.add_argument('--batch-img', type=int, default=0, help='override del batch de imagen')
    p.add_argument('--batch-txt', type=int, default=0, help='override del batch de texto')
    p.add_argument('--guardar-npz', action='store_true',
                   help='guarda tambien los .npz (por defecto NO: solo mide tiempo/RAM)')
    p.add_argument('--cpu-nombre', default=None, help='nombre del CPU a registrar (ej. "AMD Ryzen 7 5700U ...")')
    p.add_argument('--salida-dir', default='resultados_benchmark')
    return p.parse_args()


def main():
    args = parse_args()

    # --- fijar hilos ANTES de importar numpy/torch (respetan estas env vars al importar) ---
    import psutil
    hilos = args.hilos or (psutil.cpu_count(logical=False) or psutil.cpu_count(logical=True))
    for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[v] = str(hilos)

    import numpy as np
    import pyarrow       # <--- AÑADE ESTO AQUÍ
    import datasets      # <--- AÑADE ESTO AQUÍ
    import torch

    _orig_getattr = torch.nn.Module.__getattr__
    def _patched_getattr(self, name):
        if name == 'all_tied_weights_keys':
            return {}
        return _orig_getattr(self, name)
    torch.nn.Module.__getattr__ = _patched_getattr

    torch.set_num_threads(hilos)
    torch.manual_seed(123)
    np.random.seed(123)

    import bench_common as bc
    import bench_modelos as bm

    # --- resolver device y dtype ---
    device = args.device
    if device == 'auto':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    if device == 'cuda' and not torch.cuda.is_available():
        raise SystemExit('Pediste --device cuda pero torch no ve ninguna GPU.')
    en_gpu = (device == 'cuda')
    # dtype por defecto: en GPU, float16 para Qwen/InternVL (= tus notebooks de Colab, embeddings
    # identicos a los canonicos); en CPU siempre float32. CLIP/BLIP van en float32 en ambos.
    if args.dtype:
        dtype = args.dtype
    elif en_gpu and args.modelo in ('qwen_vl', 'internvl'):
        dtype = 'float16'
    else:
        dtype = 'float32'

    os.makedirs(args.salida_dir, exist_ok=True)
    host = bc.platform.node().replace(' ', '_')
    ruta_csv = os.path.join(args.salida_dir, f'benchmark_{args.modelo}_{device}_{host}.csv')
    ruta_env = os.path.join(args.salida_dir, f'entorno_{args.modelo}_{device}_{host}.json')
    dir_npz = os.path.join(args.salida_dir, 'embeddings')

    # --- entorno / reproducibilidad ---
    entorno = bc.capturar_entorno(cpu_nombre=args.cpu_nombre)
    entorno['modelo'] = args.modelo
    entorno['device'] = device
    entorno['dtype'] = dtype
    entorno['n_hilos'] = hilos
    entorno['max_n'] = args.max_n
    print('== ENTORNO ==')
    for k, v in entorno.items():
        print(f'  {k}: {v}')
    print('=============')

    # --- que datasets correr ---
    def filtra(registro):
        items = list(registro.items())
        if args.grupo != 'todos':
            items = [(n, i) for n, i in items if i['grupo'] == args.grupo]
        if args.solo:
            quiere = {s.strip() for s in args.solo.split(',') if s.strip()}
            items = [(n, i) for n, i in items if n in quiere]
        return dict(items)

    trabajos = []  # (modalidad, embed_key, es_texto, registro)
    if args.modalidad in ('imagen', 'ambas'):
        trabajos.append(('imagen', 'embed_images', False, filtra(bc.REGISTRO_IMG)))
    if args.modalidad in ('texto', 'ambas'):
        trabajos.append(('texto', 'embed_texts', True, filtra(bc.REGISTRO_TXT)))

    hechos = bc.leer_hechos(ruta_csv)
    total = sum(len(r) for _, _, _, r in trabajos)
    pendientes = sum(1 for mod, _, _, r in trabajos for n in r if (mod, n) not in hechos)
    print(f'Datasets: {total} en total, {len(hechos)} ya hechos, {pendientes} pendientes.\n')

    # --- muestreador de RAM ---
    #import datasets  # <--- AÑADE ESTA LÍNEA EXACTAMENTE AQUÍ
    #sampler = bc.RAMSampler(intervalo=0.1)
    #sampler.start()

    # --- carga del modelo (una sola vez) ---
    sampler = bc.RAMSampler(intervalo=0.1)
    print(f'Cargando modelo "{args.modelo}" en {device} (dtype={dtype})...')
    ram_antes = sampler.actual_MB()
    t0 = time.perf_counter()
    M = bm.cargar_modelo(args.modelo, device=device, dtype=dtype,
                         batch_img=args.batch_img or None, batch_txt=args.batch_txt or None)
    if en_gpu:
        torch.cuda.synchronize()
    t_carga_modelo = time.perf_counter() - t0
    ram_modelo = sampler.actual_MB() - ram_antes
    entorno.update({'model_name': M['model_name'], 'dtype': M['dtype'],
                    'batch_img': M['batch_img'], 'batch_txt': M['batch_txt'],
                    't_carga_modelo_s': round(t_carga_modelo, 2),
                    'ram_modelo_MB': round(ram_modelo, 1)})
    bc.guardar_entorno(ruta_env, entorno)
    print(f'  modelo listo en {t_carga_modelo:.1f}s | RAM del modelo ~{ram_modelo:.0f} MB '
          f'| dtype={M["dtype"]} | batch img/txt={M["batch_img"]}/{M["batch_txt"]}\n')

    hecho = 0
    for modalidad, embed_key, es_texto, registro in trabajos:
        embed_fn = M[embed_key]
        for nombre, info in registro.items():
            if (modalidad, nombre) in hechos:
                continue
            hecho += 1
            print(f'[{hecho}/{pendientes}] {modalidad}:{nombre}  (grupo {info["grupo"]}, k={info["k"]})')
            fila = {'modelo': args.modelo, 'modalidad': modalidad, 'grupo': info['grupo'],
                    'dataset': nombre, 'k_obj': info['k'], 'device': device, 'dtype': M['dtype'],
                    'batch': M['batch_img'] if not es_texto else M['batch_txt'],
                    'n_hilos': hilos, 'fecha': time.strftime('%Y-%m-%d %H:%M:%S')}
            try:
                # 1) cargar dataset (descarga la 1a vez; NO cuenta como tiempo de embedding)
                t1 = time.perf_counter()
                datos, y, names = info['loader']()
                if es_texto:
                    datos, y = bc.limpiar_textos(datos, y)
                # 2) tope opcional de N (timing representativo)
                if args.max_n and len(datos) > args.max_n:
                    idx = bc.submuestreo_estratificado(y, args.max_n, seed=123)
                    datos = [datos[i] for i in idx]; y = y[idx]
                    print(f'  submuestreo --max-n {args.max_n}: N={len(datos)}')
                t_datos = time.perf_counter() - t1
                N = len(datos)

                # 3) EMBEDDING puro (lo que reportas en la tesis)
                #    con reintento automatico si hay CUDA out of memory
                sampler = bc.RAMSampler(intervalo=0.1)
                sampler.start()
                sampler.reset_pico()
                ram_base = sampler.actual_MB()
                b_actual = fila['batch']
                while True:
                    try:
                        if en_gpu:
                            torch.cuda.synchronize(); torch.cuda.reset_peak_memory_stats()
                        t2 = time.perf_counter()
                        Xr = embed_fn(datos, batch=b_actual)
                        if en_gpu:
                            torch.cuda.synchronize()   # el GPU es asincrono: hay que esperar
                        t_embed = time.perf_counter() - t2
                        break
                    except RuntimeError as e:
                        if en_gpu and 'out of memory' in str(e).lower() and b_actual > 1:
                            torch.cuda.empty_cache()
                            b_actual = max(1, b_actual // 2)
                            print(f'    CUDA OOM -> reintento con batch={b_actual}')
                            continue
                        raise
                ram_pico = sampler.pico_MB()
                sampler.stop()
                vram_pico = (torch.cuda.max_memory_allocated() / (1024 ** 2)) if en_gpu else ''
                X = bc.normaliza_l2(Xr)
                fila['batch'] = b_actual   # el batch que finalmente funciono

                X, y = bc.limpiar_embeddings(X, y)
                k_real = int(len(np.unique(y)))
                dim = int(X.shape[1])

                if args.guardar_npz:
                    bc.guardar_npz(dir_npz, args.modelo, modalidad, nombre, X, y, names,
                                   info['grupo'], info['k'])

                fila.update({
                    'k_real': k_real, 'N': N, 'dim': dim,
                    't_datos_s': round(t_datos, 3), 't_embed_s': round(t_embed, 3),
                    'seg_por_item': round(t_embed / N, 5) if N else '',
                    'items_por_s': round(N / t_embed, 2) if t_embed else '',
                    'ram_base_MB': round(ram_base, 1), 'ram_pico_MB': round(ram_pico, 1),
                    'ram_delta_MB': round(ram_pico - ram_base, 1),
                    'vram_pico_MB': round(vram_pico, 1) if en_gpu else '',
                    'estado': 'ok', 'error': ''})
                print(f'  OK  N={N} dim={dim} | embed={t_embed:.1f}s '
                      f'({N/t_embed:.1f} it/s) | RAM pico={ram_pico:.0f} MB '
                      f'(+{ram_pico-ram_base:.0f})\n')

            except Exception as e:
                if 'sampler' in locals():
                    sampler.stop()
                fila.update({'estado': 'error', 'error': str(e)[:300]})
                print(f'  ERROR: {e}\n')

            bc.agregar_fila(ruta_csv, fila)

    print(f'\nListo. Resultados en:\n  {ruta_csv}\n  {ruta_env}')
    print('Para juntar las 3 maquinas y armar tablas:  python resumen_benchmark.py')


if __name__ == '__main__':
    main()
