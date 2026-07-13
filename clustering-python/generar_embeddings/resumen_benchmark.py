# -*- coding: utf-8 -*-
"""
resumen_benchmark.py
====================
Junta todos los CSV benchmark_*.csv (de tus 3 maquinas) y produce:

  1) Tabla de TIEMPO de embedding por dataset x modelo, una por modalidad,
     en Markdown (para revisar) y en LaTeX estilo IEEE (para pegar en la tesis).
  2) Resumen por modelo: tiempo total, throughput medio, pico de RAM.

Metrica configurable con --metrica:
  t_embed_s    -> segundos de embedding por dataset (por defecto)
  items_por_s  -> throughput (instancias/s), mas comparable entre modelos
  ram_pico_MB  -> pico de RAM

Uso:
  python resumen_benchmark.py
  python resumen_benchmark.py --metrica items_por_s
  python resumen_benchmark.py --dir resultados_benchmark --salida tablas
"""

import argparse
import glob
import os

import pandas as pd

ORDEN_MODELOS = ['clip', 'blip', 'qwen_vl', 'internvl']
NOMBRE_BONITO = {'clip': 'CLIP', 'blip': 'BLIP', 'qwen_vl': 'Qwen2-VL', 'internvl': 'InternVL2'}


def cargar(dir_):
    files = glob.glob(os.path.join(dir_, 'benchmark_*.csv'))
    if not files:
        raise SystemExit(f'No encontre benchmark_*.csv en "{dir_}".')
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df[df['estado'] == 'ok'].copy()
    if 'device' not in df.columns:      # CSV viejos sin columna device
        df['device'] = 'cpu'
    df['device'] = df['device'].fillna('cpu')
    # dedup incluyendo device: asi una misma (modelo,dataset) en CPU y en GPU coexisten
    antes = len(df)
    df = df.drop_duplicates(subset=['modelo', 'modalidad', 'dataset', 'device'], keep='last')
    if len(df) < antes:
        print(f'Aviso: {antes-len(df)} filas duplicadas (mismo dataset+device en >1 maquina); me quede con la ultima.')
    for c in ['t_embed_s', 'items_por_s', 'ram_pico_MB', 'vram_pico_MB', 'N', 'dim']:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


def pivote(df, modalidad, metrica):
    sub = df[df['modalidad'] == modalidad]
    piv = sub.pivot_table(index=['dataset'], columns='modelo', values=metrica, aggfunc='first')
    # columna N (deberia ser igual entre modelos; tomo la primera disponible)
    Ns = sub.pivot_table(index='dataset', columns='modelo', values='N', aggfunc='first')
    piv.insert(0, 'N', Ns.bfill(axis=1).iloc[:, 0].astype('Int64'))
    cols = [m for m in ORDEN_MODELOS if m in piv.columns]
    piv = piv[['N'] + cols]
    return piv.sort_values('N')


def a_markdown(piv, metrica):
    ren = {m: NOMBRE_BONITO.get(m, m) for m in piv.columns}
    p = piv.rename(columns=ren).round(2)
    return p.to_markdown()


def a_latex(piv, modalidad, metrica, etiqueta=''):
    ren = {m: NOMBRE_BONITO.get(m, m) for m in piv.columns if m in NOMBRE_BONITO}
    p = piv.rename(columns=ren)
    modelos = [c for c in p.columns if c != 'N']
    hw = {'cpu': 'en CPU', 'cuda': 'en GPU'}.get(etiqueta, '')
    esc = lambda s: str(s).replace('_', r'\_')  # los _ rompen LaTeX en modo texto
    filas = []
    for ds, row in p.iterrows():
        vals = []
        for m in modelos:
            v = row[m]
            vals.append('--' if pd.isna(v) else (f'{v:.1f}' if metrica != 'ram_pico_MB' else f'{v:.0f}'))
        n = '' if pd.isna(row['N']) else f'{int(row["N"])}'
        filas.append(f'{esc(ds)} & {n} & ' + ' & '.join(vals) + r' \\')
    col_fmt = 'l r ' + 'r ' * len(modelos)
    cap = {'t_embed_s': 'Tiempo de generacion de embeddings (s)',
           'items_por_s': 'Throughput de generacion de embeddings (instancias/s)',
           'ram_pico_MB': 'Pico de RAM (MB) durante la generacion de embeddings',
           'vram_pico_MB': 'Pico de VRAM (MB) durante la generacion de embeddings'}.get(metrica, metrica)
    cap = f'{cap} {hw} ({modalidad}).'.replace('  ', ' ')
    lab = f'tab:{metrica}_{modalidad}' + (f'_{etiqueta}' if etiqueta else '')
    head = ' & '.join([r'\textbf{Dataset}', r'\textbf{N}'] + [rf'\textbf{{{m}}}' for m in modelos])
    return (
        '\\begin{table}[t]\n\\centering\n'
        f'\\caption{{{cap}}}\n'
        f'\\label{{{lab}}}\n'
        f'\\begin{{tabular}}{{{col_fmt.strip()}}}\n\\toprule\n'
        f'{head} \\\\\n\\midrule\n'
        + '\n'.join(filas) +
        '\n\\bottomrule\n\\end{tabular}\n\\end{table}\n'
    )


def resumen_modelos(df):
    g = df.groupby('modelo')
    d = {
        'datasets': g.size(),
        't_embed_total_s': g['t_embed_s'].sum().round(1),
        'it_por_s_medio': g['items_por_s'].mean().round(1),
        'ram_pico_max_MB': g['ram_pico_MB'].max().round(0),
    }
    if 'vram_pico_MB' in df.columns and df['vram_pico_MB'].notna().any():
        d['vram_pico_max_MB'] = g['vram_pico_MB'].max().round(0)
    r = pd.DataFrame(d)
    return r.reindex([m for m in ORDEN_MODELOS if m in r.index])


def procesar(df, salida, metrica, etiqueta):
    """Genera resumen + tablas por modalidad para un subconjunto (un device)."""
    suf = f'_{etiqueta}' if etiqueta else ''
    print(f'\n===== RESUMEN POR MODELO{(" ["+etiqueta+"]") if etiqueta else ""} =====')
    rm = resumen_modelos(df)
    print(rm.to_markdown())
    rm.to_csv(os.path.join(salida, f'resumen_por_modelo{suf}.csv'))

    for modalidad in ['imagen', 'texto']:
        if not (df['modalidad'] == modalidad).any():
            continue
        piv = pivote(df, modalidad, metrica)
        print(f'\n===== {modalidad.upper()} — {metrica}{(" ["+etiqueta+"]") if etiqueta else ""} =====')
        print(a_markdown(piv, metrica))
        tex = a_latex(piv, modalidad, metrica, etiqueta)
        base = f'tabla_{metrica}_{modalidad}{suf}'
        with open(os.path.join(salida, base + '.tex'), 'w', encoding='utf-8') as f:
            f.write(tex)
        piv.to_csv(os.path.join(salida, base + '.csv'))
        print(f'  -> {os.path.join(salida, base + ".tex")}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='resultados_benchmark')
    ap.add_argument('--salida', default='tablas')
    ap.add_argument('--metrica', default='t_embed_s',
                    choices=['t_embed_s', 'items_por_s', 'ram_pico_MB', 'vram_pico_MB'])
    ap.add_argument('--device', default='all', choices=['all', 'cpu', 'cuda'],
                    help='filtra por device; "all" separa las tablas por cada device presente')
    args = ap.parse_args()

    df = cargar(args.dir)
    os.makedirs(args.salida, exist_ok=True)

    if args.device != 'all':
        df = df[df['device'] == args.device]
        procesar(df, args.salida, args.metrica, args.device)
    else:
        devices = sorted(df['device'].dropna().unique())
        if len(devices) <= 1:
            procesar(df, args.salida, args.metrica, devices[0] if devices else '')
        else:
            print(f'Detecte {len(devices)} devices: {devices}. Genero tablas separadas por cada uno.')
            for d in devices:
                procesar(df[df['device'] == d], args.salida, args.metrica, d)

    print('\nLaTeX generado requiere \\usepackage{booktabs}.')


if __name__ == '__main__':
    main()
