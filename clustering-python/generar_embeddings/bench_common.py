# -*- coding: utf-8 -*-
"""
bench_common.py
================
Utilidades compartidas para el benchmark de generacion de embeddings en CPU.

Contiene, VERBATIM desde tus notebooks de Colab (para que N, dimension y
cardinalidad coincidan exactamente con lo ya generado en GPU):
  - REGISTRO_IMG (17 datasets) y REGISTRO_TXT (15 datasets)
  - todos los loaders de imagen y texto
  - limpiar_textos / limpiar_embeddings / normaliza_l2

Y agrega lo especifico del benchmark en CPU:
  - RAMSampler       : mide el pico de RAM (RSS) del proceso durante una ventana
  - capturar_entorno : huella de hardware/software para la seccion de reproducibilidad
  - submuestreo_estratificado : para el modo --max-n (timing representativo sin correr dias)
  - helpers de CSV y guardado opcional de .npz

Dependencias: numpy, Pillow, psutil, y (perezosas, segun dataset) scikit-learn,
datasets, medmnist, nltk, torchvision.
"""

import os
import time
import csv
import json
import platform
import threading
import numpy as np
from PIL import Image

# ---------------------------------------------------------------------------
# 1) UTILIDADES DE LIMPIEZA / NORMALIZACION  (identicas a tus notebooks)
# ---------------------------------------------------------------------------
def limpiar_textos(textos, labels, min_chars=10):
    ok_t, ok_y, vistos = [], [], set()
    for t, l in zip(textos, labels):
        tn = ' '.join(str(t).split())
        if len(tn) < min_chars or tn in vistos:
            continue
        vistos.add(tn); ok_t.append(t); ok_y.append(l)
    print(f'  limpiar_textos: {len(textos)} -> {len(ok_t)} (-{len(textos)-len(ok_t)})')
    return ok_t, np.array(ok_y)


def limpiar_embeddings(X, y, min_norm=1e-8):
    n0 = len(X)
    keep = np.linalg.norm(X, axis=1) > min_norm
    X, y = X[keep], y[keep]
    _, idx = np.unique(X, axis=0, return_index=True)
    idx = np.sort(idx)
    X, y = X[idx], y[idx]
    print(f'  limpiar_embeddings: {n0} -> {len(X)} (norma0/duplicados: -{n0-len(X)})')
    return X, y


def normaliza_l2(X):
    n = np.linalg.norm(X, axis=1, keepdims=True)
    return (X / np.maximum(n, 1e-10)).astype(np.float32)


def guardar_npz(base_dir, modelo, modalidad, nombre, X, y, names, grupo, k):
    """Guarda el .npz en el MISMO formato que tus notebooks (opcional en el benchmark)."""
    if names is None:
        names = [str(c) for c in sorted(np.unique(y))]
    out = os.path.join(base_dir, modelo, modalidad)
    os.makedirs(out, exist_ok=True)
    np.savez_compressed(
        os.path.join(out, nombre + '.npz'),
        X=X, y=y.astype(int),
        class_names=np.array(names, dtype=object),
        grupo=grupo, k=k, modalidad=modalidad, modelo=modelo)


# ---------------------------------------------------------------------------
# 2) LOADERS DE IMAGEN  -> (lista_PIL, y(int), class_names|None)
# ---------------------------------------------------------------------------
def _to_pil_rgb(arr):
    a = np.asarray(arr)
    if a.dtype != np.uint8:
        a = a.astype(np.float64); a = a - a.min()
        m = a.max() if a.max() > 0 else 1.0
        a = (255.0 * a / m).astype(np.uint8)
    return Image.fromarray(a).convert('RGB')


def load_olivetti():
    from sklearn.datasets import fetch_olivetti_faces
    d = fetch_olivetti_faces()
    return [_to_pil_rgb(i) for i in d.images], d.target.astype(int), [f'p{i}' for i in range(40)]


def load_digits_sklearn():
    from sklearn.datasets import load_digits
    d = load_digits()
    return [_to_pil_rgb(i) for i in d.images], d.target.astype(int), [str(i) for i in range(10)]


def load_lfw():
    from sklearn.datasets import fetch_lfw_people
    d = fetch_lfw_people(min_faces_per_person=30, resize=0.5, color=False)
    return [_to_pil_rgb(i) for i in d.images], d.target.astype(int), list(d.target_names)


def _tv(cls, **kw):
    import torchvision
    ds = getattr(torchvision.datasets, cls)(root='data', download=True, **kw)
    im, y = [], []
    for a, b in ds:
        im.append(a.convert('RGB')); y.append(int(b))
    return im, np.array(y)


def load_semeion():      i, y = _tv('SEMEION');                 return i, y, [str(j) for j in range(10)]
def load_mnist_test():   i, y = _tv('MNIST', train=False);      return i, y, [str(j) for j in range(10)]
def load_fashion_test(): i, y = _tv('FashionMNIST', train=False); return i, y, ['Tshirt','Trouser','Pullover','Dress','Coat','Sandal','Shirt','Sneaker','Bag','Boot']
def load_cifar10_test(): i, y = _tv('CIFAR10', train=False);    return i, y, ['airplane','auto','bird','cat','deer','dog','frog','horse','ship','truck']
def load_stl10_train():  i, y = _tv('STL10', split='train');    return i, y, ['airplane','bird','car','cat','deer','dog','horse','monkey','ship','truck']


def load_usps_full():
    ai, ay = _tv('USPS', train=True); bi, by = _tv('USPS', train=False)
    return ai + bi, np.concatenate([ay, by]), [str(j) for j in range(10)]


def load_pet_full():
    ai, ay = _tv('OxfordIIITPet', split='trainval', target_types='category')
    bi, by = _tv('OxfordIIITPet', split='test', target_types='category')
    return ai + bi, np.concatenate([ay, by]), [f'breed{j}' for j in range(37)]


CIFAR100_COARSE = np.array([4,1,14,8,0,6,7,7,18,3,3,14,9,18,7,11,3,9,7,11,6,11,5,10,7,6,13,15,3,15,0,11,1,10,12,14,16,9,11,5,5,19,8,8,15,13,14,17,18,10,16,4,17,4,2,0,17,4,18,17,10,3,2,12,12,16,12,1,9,19,2,10,0,1,16,12,9,13,15,13,16,19,2,4,6,19,5,5,8,19,18,1,2,15,6,0,17,8,14,13])


def load_cifar100_super_test():
    import torchvision
    ds = torchvision.datasets.CIFAR100(root='data', train=False, download=True)
    im, y = [], []
    for a, fine in ds:
        im.append(a.convert('RGB')); y.append(int(CIFAR100_COARSE[fine]))
    return im, np.array(y), [f'super{j}' for j in range(20)]


def _med(cls, splits, k):
    import medmnist
    C = getattr(medmnist, cls); im, ys = [], []
    for sp in splits:
        ds = C(split=sp, download=True)
        X = ds.imgs; Y = ds.labels.reshape(-1).astype(int)
        for j in range(len(X)):
            im.append(_to_pil_rgb(X[j])); ys.append(int(Y[j]))
    return im, np.array(ys), [str(j) for j in range(k)]


def load_retina(): return _med('RetinaMNIST', ['train', 'val', 'test'], 5)
def load_derma():  return _med('DermaMNIST', ['test'], 7)
def load_blood():  return _med('BloodMNIST', ['test'], 8)
def load_path():   return _med('PathMNIST', ['test'], 9)
def load_organc(): return _med('OrganCMNIST', ['test'], 11)
def load_organs(): return _med('OrganSMNIST', ['test'], 11)


# ---------------------------------------------------------------------------
# 3) LOADERS DE TEXTO  -> (lista_textos, y(int), class_names|None)
# ---------------------------------------------------------------------------
def _ld(*a, **k):
    from datasets import load_dataset
    try:
        return load_dataset(*a, **k)
    except Exception as e:
        msg = str(e).lower()
        if 'script' in msg or 'invalid hf uri' in msg or 'no longer supported' in msg:
            try:
                return load_dataset(*a, revision='refs/convert/parquet', **k)
            except Exception:
                raise e
        raise


def _all(d):
    from datasets import concatenate_datasets
    return concatenate_datasets([d[s] for s in d.keys()])


def _cap(textos, y, n=10000, seed=123):
    """Acota un split grande a n, ESTRATIFICADO y reproducible (seed=123)."""
    if len(textos) <= n:
        return textos, np.asarray(y)
    y = np.asarray(y); rng = np.random.RandomState(seed); idx = []
    for c in np.unique(y):
        ci = np.where(y == c)[0]
        take = max(1, int(round(n * len(ci) / len(y))))
        idx.extend(rng.choice(ci, size=min(take, len(ci)), replace=False))
    idx = np.array(sorted(idx))[:n]
    return [textos[i] for i in idx], y[idx]


def load_bbc_news():    f = _all(_ld('SetFit/bbc-news')); return list(f['text']), np.array(f['label']), None
def load_emotion_test():s = _ld('dair-ai/emotion', 'split')['test']; return list(s['text']), np.array(s['label']), None


def load_liar():
    from datasets import concatenate_datasets
    d = _ld('chengxuphd/liar2')
    f = concatenate_datasets([d['train'], d['validation'], d['test']])
    t, y = _cap(list(f['statement']), np.array(f['label']))
    return t, y, None


def load_massive_scenario(): s = _ld('AmazonScience/massive', 'en-US')['test']; return list(s['utt']), np.array(s['scenario']), None


def load_goemotions_single():
    s = _ld('google-research-datasets/go_emotions', 'simplified')['test']; t, y = [], []
    for ex in s:
        if len(ex['labels']) == 1:
            t.append(ex['text']); y.append(ex['labels'][0])
    return t, np.array(y), None


def load_tweeteval_emoji_val(): s = _ld('cardiffnlp/tweet_eval', 'emoji')['validation']; return list(s['text']), np.array(s['label']), None


def load_trec():
    s = _all(_ld('CogComp/trec'))
    lf = 'coarse_label' if 'coarse_label' in s.column_names else 'label-coarse'
    return list(s['text']), np.array(s[lf]), None


def load_sst5_train(): s = _ld('SetFit/sst5')['train']; return list(s['text']), np.array(s['label']), None


def load_20ng_test():
    from sklearn.datasets import fetch_20newsgroups
    d = fetch_20newsgroups(subset='test', remove=('headers', 'footers', 'quotes'))
    return list(d.data), np.array(d.target), list(d.target_names)


def load_reuters_r10():
    import nltk; nltk.download('reuters', quiet=True)
    from nltk.corpus import reuters
    from collections import Counter
    dc = {fid: reuters.categories(fid)[0] for fid in reuters.fileids() if len(reuters.categories(fid)) == 1}
    top = [c for c, _ in Counter(dc.values()).most_common(10)]
    n2i = {c: i for i, c in enumerate(top)}; t, y = [], []
    for fid, c in dc.items():
        if c in n2i:
            t.append(reuters.raw(fid)); y.append(n2i[c])
    return t, np.array(y), top


def load_patent(): s = _ld('ccdv/patent-classification', 'abstract')['test']; return list(s['text']), np.array(s['label']), None
def load_scotus(): s = _ld('coastalcph/lex_glue', 'scotus')['test']; return list(s['text']), np.array(s['label']), None
def load_mtop_domain(): s = _ld('mteb/mtop_domain', 'en')['test']; return list(s['text']), np.array(s['label']), None
def load_ag_news(): s = _ld('fancyzhx/ag_news')['test']; return list(s['text']), np.array(s['label']), None


def load_meld_e():
    s = _ld('boltuix/emotions-dataset')['train']
    from collections import Counter
    conteo = Counter(s['Label']); top7 = [c for c, _ in conteo.most_common(7)]
    n2i = {c: i for i, c in enumerate(top7)}; t, y = [], []
    for txt, lab in zip(s['Sentence'], s['Label']):
        if lab in n2i:
            t.append(txt); y.append(n2i[lab])
    t, y = _cap(t, np.array(y))
    return t, y, top7


def load_dbpedia14():
    s = _ld('fancyzhx/dbpedia_14')['test']
    t, y = _cap(list(s['content']), np.array(s['label']))
    return t, y, None


def load_yahoo():
    s = _ld('community-datasets/yahoo_answers_topics')['test']
    t = [f"{a} {b}" for a, b in zip(s['question_title'], s['best_answer'])]
    t, y = _cap(t, np.array(s['topic']))
    return t, y, None


# ---------------------------------------------------------------------------
# 4) REGISTROS  (identicos a tus notebooks)
# ---------------------------------------------------------------------------
REGISTRO_IMG = {
 'orl':         dict(grupo='A', k=40, loader=load_olivetti),
 'retinamnist': dict(grupo='A', k=5,  loader=load_retina),
 'digits':      dict(grupo='A', k=10, loader=load_digits_sklearn),
 'lfw':         dict(grupo='A', k=34, loader=load_lfw),
 'bloodmnist':  dict(grupo='A', k=8,  loader=load_blood),
 'pathmnist':   dict(grupo='A', k=9,  loader=load_path),
 'organcmnist': dict(grupo='A', k=11, loader=load_organc),
 'usps':        dict(grupo='A', k=10, loader=load_usps_full),
 'semeion':     dict(grupo='B', k=10, loader=load_semeion),
 'dermamnist':  dict(grupo='B', k=7,  loader=load_derma),
 'stl10':       dict(grupo='B', k=10, loader=load_stl10_train),
 'oxfordpet':   dict(grupo='B', k=37, loader=load_pet_full),
 'organsmnist': dict(grupo='B', k=11, loader=load_organs),
 'mnist':       dict(grupo='B', k=10, loader=load_mnist_test),
 'fashionmnist':dict(grupo='B', k=10, loader=load_fashion_test),
 'cifar10':     dict(grupo='B', k=10, loader=load_cifar10_test),
 'cifar100':    dict(grupo='B', k=20, loader=load_cifar100_super_test),
}

REGISTRO_TXT = {
 'ag_news':        dict(grupo='A', k=4,  loader=load_ag_news),
 'sst5':           dict(grupo='A', k=5,  loader=load_sst5_train),
 'liar':           dict(grupo='A', k=6,  loader=load_liar),
 'meld_e':         dict(grupo='A', k=7,  loader=load_meld_e),
 'yahoo':          dict(grupo='A', k=10, loader=load_yahoo),
 'dbpedia14':      dict(grupo='A', k=14, loader=load_dbpedia14),
 'goemotions':     dict(grupo='A', k=28, loader=load_goemotions_single),
 'bbc_news':       dict(grupo='B', k=5,  loader=load_bbc_news),
 'emotion':        dict(grupo='B', k=6,  loader=load_emotion_test),
 'trec':           dict(grupo='B', k=6,  loader=load_trec),
 'patent':         dict(grupo='B', k=9,  loader=load_patent),
 'reuters_r10':    dict(grupo='B', k=10, loader=load_reuters_r10),
 'scotus':         dict(grupo='B', k=13, loader=load_scotus),
 'tweeteval_emoji':dict(grupo='B', k=20, loader=load_tweeteval_emoji_val),
 '20ng':           dict(grupo='B', k=20, loader=load_20ng_test),
}


# ---------------------------------------------------------------------------
# 5) BENCHMARK: muestreador de RAM, entorno, submuestreo, CSV
# ---------------------------------------------------------------------------
class RAMSampler:
    """Muestrea el RSS del proceso en un hilo de fondo y guarda el pico.

    Uso:
        s = RAMSampler(); s.start()
        s.reset_pico()          # justo antes de la ventana que quieres medir
        ... trabajo ...
        pico = s.pico_MB()      # pico de RSS durante la ventana
        s.stop()
    """
    def __init__(self, intervalo=0.1):
        import psutil
        self._proc = psutil.Process(os.getpid())
        self._intervalo = intervalo
        self._pico = 0.0
        self._corriendo = False
        self._hilo = None

    def _rss_mb(self):
        return self._proc.memory_info().rss / (1024 ** 2)

    def _loop(self):
        while self._corriendo:
            r = self._rss_mb()
            if r > self._pico:
                self._pico = r
            time.sleep(self._intervalo)

    def start(self):
        self._corriendo = True
        self._hilo = threading.Thread(target=self._loop, daemon=True)
        self._hilo.start()

    def reset_pico(self):
        self._pico = self._rss_mb()

    def pico_MB(self):
        return max(self._pico, self._rss_mb())

    def actual_MB(self):
        return self._rss_mb()

    def stop(self):
        self._corriendo = False
        if self._hilo is not None:
            self._hilo.join(timeout=1.0)


def submuestreo_estratificado(y, n, seed=123):
    """Devuelve indices ORDENADOS para quedarse con ~n items respetando proporciones
    de clase (reproducible, seed=123). Sirve para el modo --max-n (timing rapido)."""
    y = np.asarray(y)
    if len(y) <= n:
        return np.arange(len(y))
    rng = np.random.RandomState(seed); idx = []
    for c in np.unique(y):
        ci = np.where(y == c)[0]
        take = max(1, int(round(n * len(ci) / len(y))))
        idx.extend(rng.choice(ci, size=min(take, len(ci)), replace=False))
    return np.array(sorted(idx))[:n]


def capturar_entorno(cpu_nombre=None):
    """Huella de hardware/software para la seccion de reproducibilidad de la tesis."""
    import psutil
    info = {
        'fecha': time.strftime('%Y-%m-%d %H:%M:%S'),
        'host': platform.node(),
        'so': f'{platform.system()} {platform.release()}',
        'cpu': cpu_nombre or platform.processor() or 'desconocido',
        'cpu_cores_fisicos': psutil.cpu_count(logical=False),
        'cpu_hilos_logicos': psutil.cpu_count(logical=True),
        'ram_total_GB': round(psutil.virtual_memory().total / (1024 ** 3), 1),
        'python': platform.python_version(),
        'numpy': np.__version__,
    }
    try:
        import torch; info['torch'] = torch.__version__
        info['torch_num_threads'] = torch.get_num_threads()
        info['cuda_disponible'] = bool(torch.cuda.is_available())
        info['gpu'] = torch.cuda.get_device_name(0) if torch.cuda.is_available() else None
    except Exception:
        info['torch'] = None; info['torch_num_threads'] = None
        info['cuda_disponible'] = None; info['gpu'] = None
    try:
        import transformers; info['transformers'] = transformers.__version__
    except Exception:
        info['transformers'] = None
    return info


# --- CSV ---------------------------------------------------------------------
COLUMNAS = ['modelo', 'modalidad', 'grupo', 'dataset', 'k_obj', 'k_real', 'N', 'dim',
            'device', 't_datos_s', 't_embed_s', 'seg_por_item', 'items_por_s',
            'ram_base_MB', 'ram_pico_MB', 'ram_delta_MB', 'vram_pico_MB',
            'dtype', 'batch', 'n_hilos', 'estado', 'error', 'fecha']


def leer_hechos(ruta_csv):
    """Devuelve set de (modalidad, dataset) ya calculados OK (para reanudar)."""
    hechos = set()
    if not os.path.exists(ruta_csv):
        return hechos
    with open(ruta_csv, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            if row.get('estado') == 'ok':
                hechos.add((row['modalidad'], row['dataset']))
    return hechos


def agregar_fila(ruta_csv, fila):
    """Agrega una fila (dict) al CSV, creando cabecera si hace falta."""
    nuevo = not os.path.exists(ruta_csv)
    with open(ruta_csv, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=COLUMNAS)
        if nuevo:
            w.writeheader()
        w.writerow({c: fila.get(c, '') for c in COLUMNAS})


def guardar_entorno(ruta_json, info):
    with open(ruta_json, 'w', encoding='utf-8') as f:
        json.dump(info, f, ensure_ascii=False, indent=2)
