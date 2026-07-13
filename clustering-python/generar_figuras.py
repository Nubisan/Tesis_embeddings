# -*- coding: utf-8 -*-
"""
Genera TODAS las figuras de la seccion de Resultados y Discusion a partir de:
  - resultados_largo.csv  (salida de build_results_table.py)
  - peakRAM_log.csv       (log de memoria)
  - runner_log.csv        (log de tiempo reloj)

Produce, en la carpeta SALIDA, cada figura en PNG:
  1. crossmodal_externa          -> comparacion de modelos de embedding x modalidad
  2. boxplot_<MET>_<modalidad>   -> distribucion por algoritmo, facetado por modelo
  3. victorias_externa           -> conteo de mejores resultados (ARI/AMI/NMI)
  4. boxplot_silhouette_<mod>    -> distribucion de Silhouette
  5. victorias_silhouette        -> conteo de maximos en validacion interna
  6. tiempo                      -> tiempo medio de ejecucion por algoritmo
  7. ram                         -> pico de RAM por algoritmo
  8. tiempo vs n                 -> comparacion del tiempo con instancias
  9. tiempo reloj               -> comparacion del tiempo reloj x modelos de embedding

Ejecutar:    python generar_figuras.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ============================================================
# CONFIGURACION
# ============================================================
ARCHIVO_RESULTADOS = "resultados_largo.csv"
ARCHIVO_RAM        = "peakRAM_log.csv"
ARCHIVO_RUNNER     = "runner_log.csv"   # tiempo reloj por dataset
SALIDA             = "figuras"          # carpeta de salida
DPI                = 300

# Layout de los boxplots facetados por modelo:
#   "1x4" -> 4 paneles en fila    (para figura a dos columnas: \begin{figure*})
#   "2x2" -> 4 paneles en rejilla (para figura de una sola columna)
LAYOUT = "1x4"

FUENTE = "serif"
rcParams["font.family"] = "serif" if FUENTE == "serif" else "sans-serif"
if FUENTE != "serif":
    rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
rcParams["font.size"] = 10

rcParams["axes.spines.top"]   = False
rcParams["axes.spines.right"] = False

# Orden de algoritmos
ORDEN_ALGOS = ["CSCLP", "MILP-KM", "SCK", "HCAKC",
               "K-MeansBA", "K-MeansACO", "K-MeansPSO", "K-MedoidsSC"]
COLOR_ALGO = {
    "CSCLP":       "#0072B2",  # azul
    "MILP-KM":     "#E69F00",  # naranja
    "SCK":         "#009E73",  # verde
    "HCAKC":       "#CC79A7",  # rosa
    "K-MeansBA":   "#56B4E9",  # celeste
    "K-MeansACO":  "#D55E00",  # bermellon
    "K-MeansPSO":  "#F0E442",  # amarillo
    "K-MedoidsSC": "#785EF0",  # violeta
}
ORDEN_MODELOS = ["CLIP", "BLIP", "Qwen2-VL", "InternVL2"]
# Colores por modalidad, coherentes con victorias externas (azul / naranja)
COL_IMG, COL_TXT = "#0072B2", "#E69F00"
# Borde gris oscuro fino en cajas y marcadores
EDGE = dict(edgecolor="#333333", linewidth=0.6)

# Algoritmos que se excluyen del promedio "competitivo" del cross-modal
NO_COMPETITIVOS = ["K-MeansPSO", "K-MedoidsSC"]


# ================
# CARGA DE DATOS
# ================
def cargar():
    r = pd.read_csv(ARCHIVO_RESULTADOS)
    ram = None
    if os.path.exists(ARCHIVO_RAM):
        ram = pd.read_csv(ARCHIVO_RAM)
        ram["modelo"] = ram["model"].map(
            {"clip": "CLIP", "blip": "BLIP", "internvl": "InternVL2", "qwen_vl": "Qwen2-VL"})
        ram["algoritmo"] = ram["algorithm"].map(
            {"ACO": "K-MeansACO", "BAT": "K-MeansBA", "CSCLP": "CSCLP", "HCAKC": "HCAKC",
             "KM-MILP": "MILP-KM", "Kmedoids": "K-MedoidsSC", "PSO": "K-MeansPSO", "SCK1": "SCK"})
    return r, ram


def guardar(fig, nombre):
    fig.savefig(os.path.join(SALIDA, nombre + ".png"), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("  ->", nombre)


# =====================================
# 1) CROSS-MODAL (validacion externa)
# =====================================
def fig_crossmodal(r):
    rc = r[~r.algoritmo.isin(NO_COMPETITIVOS)]
    mets = ["ARI", "AMI", "NMI"]
    modelos = ORDEN_MODELOS
    y = np.arange(len(modelos))[::-1]
    fig, axes = plt.subplots(3, 1, figsize=(6, 8))
    for ax, met in zip(axes, mets):
        p = rc.pivot_table(index="modelo", columns="modalidad",
                           values=met, aggfunc="mean").reindex(modelos)
        for yi, m in zip(y, modelos):
            ax.plot([p.loc[m, "imagen"], p.loc[m, "texto"]], [yi, yi],
                    color="#BBBBBB", lw=2, zorder=1)
        ax.scatter(p["imagen"], y, s=70, color=COL_IMG, label="Imagen",
                   zorder=3, **EDGE)
        ax.scatter(p["texto"], y, s=70, color=COL_TXT, label="Texto",
                   zorder=3, **EDGE)
        ax.set_title(met, fontsize=11)
        ax.set_yticks(y); ax.set_yticklabels(modelos, fontsize=9)
        ax.grid(axis="x", ls=":", alpha=0.5); ax.set_axisbelow(True)
        ax.set_xlabel("Valor promedio", fontsize=9)
        ax.margins(y=0.15)
    axes[1].legend(loc="upper right", frameon=True, fontsize=8.5)
    fig.tight_layout()
    guardar(fig, "crossmodal_externa")

# ============================================================
# 2) BOXPLOTS por metrica y modalidad, facetado por modelo
# ============================================================
def fig_boxplot(r, metrica, modalidad, nombre):
    sub = r[r.modalidad == modalidad]
    if LAYOUT == "2x2":
        fig, axes = plt.subplots(2, 2, figsize=(6.6, 6.2), sharey=True)
        axes = axes.ravel()
        rot, fs = 90, 7.0
    else:  # "1x4"
        fig, axes = plt.subplots(1, 4, figsize=(13, 3.8), sharey=True)
        rot, fs = 90, 7.5
    for ax, mod in zip(axes, ORDEN_MODELOS):
        s = sub[sub.modelo == mod]
        dat = [s[s.algoritmo == a][metrica].dropna().values for a in ORDEN_ALGOS]
        bp = ax.boxplot(dat, patch_artist=True, widths=0.6)
        for patch, a in zip(bp["boxes"], ORDEN_ALGOS):
            patch.set_facecolor(COLOR_ALGO[a]); patch.set_alpha(0.85)
            patch.set_edgecolor("#333333"); patch.set_linewidth(0.6)
        for med in bp["medians"]:
            med.set_color("#333333"); med.set_linewidth(1.3)
        ax.set_title(mod, fontsize=11)
        ax.set_xticks(range(1, len(ORDEN_ALGOS) + 1))
        ax.set_xticklabels(ORDEN_ALGOS, rotation=rot, fontsize=fs)
        ax.axhline(0, color="gray", ls="--", lw=0.7)
        ax.grid(axis="y", ls=":", alpha=0.5); ax.set_axisbelow(True)
    # etiqueta del eje Y: en 1x4 solo el primero; en 2x2 los de la izquierda
    etiqueta = metrica.replace("_mean", "")
    if LAYOUT == "2x2":
        axes[0].set_ylabel(etiqueta); axes[2].set_ylabel(etiqueta)
    else:
        axes[0].set_ylabel(etiqueta)
    fig.tight_layout()
    guardar(fig, nombre)


# ============================================================
# 3 y 5) CONTEO DE VICTORIAS (mejor algoritmo por dataset x modelo)
# ============================================================
def conteo_victorias(r, metricas):
    """Para cada (dataset, modelo) cuenta que algoritmo tiene el maximo en cada metrica."""
    out = {m: {a: 0 for a in ORDEN_ALGOS} for m in metricas}
    for (ds, mo), g in r.groupby(["name", "modelo"]):
        for m in metricas:
            gg = g.dropna(subset=[m])
            if len(gg):
                gan = gg.loc[gg[m].idxmax(), "algoritmo"]
                if gan in out[m]:
                    out[m][gan] += 1
    return out


def fig_victorias_externa(r):
    metricas = ["ARI", "AMI", "NMI"]
    out = conteo_victorias(r, metricas)
    x = np.arange(len(ORDEN_ALGOS)); w = 0.26
    fig, ax = plt.subplots(figsize=(9, 3.8))
    cols = {"ARI": "#0072B2", "AMI": "#E69F00", "NMI": "#009E73"}
    for i, m in enumerate(metricas):
        vals = [out[m][a] for a in ORDEN_ALGOS]
        ax.bar(x + (i - 1) * w, vals, w, label=m, color=cols[m], **EDGE)
    ax.set_xticks(x); ax.set_xticklabels(ORDEN_ALGOS, rotation=30, ha="right", fontsize=8.5)
    ax.set_ylabel("Frecuencia"); ax.legend(frameon=True)
    ax.grid(axis="y", ls=":", alpha=0.5); ax.set_axisbelow(True)
    fig.tight_layout()
    guardar(fig, "victorias_externa")


def fig_victorias_silhouette(r):
    out = conteo_victorias(r, ["Silhouette_mean"])["Silhouette_mean"]
    orden = sorted(ORDEN_ALGOS, key=lambda a: out[a], reverse=True)
    vals = [out[a] for a in orden]
    fig, ax = plt.subplots(figsize=(8, 3.8))
    ax.bar(range(len(orden)), vals, color=[COLOR_ALGO[a] for a in orden], alpha=0.85, **EDGE)
    ax.set_xticks(range(len(orden))); ax.set_xticklabels(orden, rotation=30, ha="right", fontsize=8.5)
    ax.set_ylabel("Frecuencia S(i)")
    ax.grid(axis="y", ls=":", alpha=0.5); ax.set_axisbelow(True)
    fig.tight_layout()
    guardar(fig, "victorias_silhouette")


# ========================
# 6) TIEMPO DE EJECUCION
# ========================
def fig_tiempo(r):
    from matplotlib import ticker
    med = r.groupby("algoritmo")["Execution_Time"].mean().reindex(ORDEN_ALGOS)
    orden = med.sort_values(ascending=True)
    y = np.arange(len(orden))
    fig, ax = plt.subplots(figsize=(8, 4))
    for yi, (a, v) in zip(y, orden.items()):
        ax.plot([0, v], [yi, yi], color="#DDDDDD", lw=1, zorder=1)
    ax.scatter(orden.values, y, s=90,
               color=[COLOR_ALGO[a] for a in orden.index], zorder=3, **EDGE)
    ax.set_yticks(y); ax.set_yticklabels(orden.index, fontsize=9)
    ax.set_xlabel("Tiempo promedio de ejecución (s)")
    ax.set_xlim(left=0)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x)}"))
    ax.grid(axis="x", ls=":", alpha=0.5); ax.set_axisbelow(True)
    ax.margins(y=0.08)
    fig.tight_layout()
    guardar(fig, "tiempo")

# ====================================================
# 7) RAM
#    AGG_RAM: "max" (pico peor caso), "mean" o "sum".
# ====================================================
AGG_RAM = "max"


def fig_ram(ram):
    if ram is None:
        print("  [aviso] no hay log de RAM, se omite figura de RAM")
        return
    p = (ram.groupby(["algoritmo", "modalidad"])["peak_ram_MB"]
            .agg(AGG_RAM).unstack().reindex(ORDEN_ALGOS))
    x = np.arange(len(ORDEN_ALGOS))
    etiqueta = {"max": "Pico máximo de RAM (MB)",
                "mean": "Pico promedio de RAM (MB)",
                "sum": "RAM total (MB)"}.get(AGG_RAM, "RAM (MB)")
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(x, p.get("imagen"), marker="o", ms=7, lw=1.8, color=COL_IMG,
            label="Imagen", mec="#333333", mew=0.5)
    ax.plot(x, p.get("texto"), marker="s", ms=7, lw=1.8, color=COL_TXT,
            label="Texto", mec="#333333", mew=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(ORDEN_ALGOS, rotation=30, ha="right", fontsize=8.5)
    ax.set_ylabel(etiqueta)
    ax.set_ylim(bottom=0)
    ax.grid(ls=":", alpha=0.5); ax.set_axisbelow(True)
    ax.legend(frameon=True)
    fig.tight_layout()
    guardar(fig, "ram")

# ==================================
# 8) TIEMPO vs NUMERO DE INSTANCIAS
# ==================================
def fig_tiempo_por_tamano(r):
    bins = [0, 2000, 7000, np.inf]
    labels = ["Pequeños\n(<2000)", "Medianos\n(2000–7000)", "Grandes\n(>7000)"]
    r = r.copy()
    r["grupo"] = pd.cut(r["n"], bins=bins, labels=labels)
    g = r.groupby("grupo", observed=True)["Execution_Time"].mean()
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.bar(range(len(g)), g.values, width=0.6,
           color=["#56B4E9", "#0072B2", "#D55E00"], **EDGE)
    ax.set_xticks(range(len(g))); ax.set_xticklabels(g.index, fontsize=9)
    ax.set_ylabel("Tiempo promedio de clustering (s)")
    ax.set_xlabel("Tamaño del conjunto de datos (n)")
    ax.grid(axis="y", ls=":", alpha=0.5); ax.set_axisbelow(True)
    for i, v in enumerate(g.values):
        ax.text(i, v + 5, f"{v:.0f} s", ha="center", fontsize=9)
    fig.tight_layout()
    guardar(fig, "tiempo_por_tamano")


# =========================================
# 9) TIEMPO RELOJ POR MODELO Y MODALIDAD
# =========================================
def fig_tiempo_por_modelo():
    if not os.path.exists(ARCHIVO_RUNNER):
        print("  [aviso] no hay runner_log, se omite tiempo_por_modelo")
        return
    rl = pd.read_csv(ARCHIVO_RUNNER)
    rl["modelo"] = rl["model"].map({"clip": "CLIP", "blip": "BLIP",
                                    "qwen_vl": "Qwen2-VL", "internvl": "InternVL2"})
    t = (rl.groupby(["modelo", "modalidad"])["elapsed_s"].sum() / 3600).unstack()
    t = t.reindex(ORDEN_MODELOS)
    y = np.arange(len(ORDEN_MODELOS))[::-1]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    for yi, m in zip(y, ORDEN_MODELOS):
        ax.plot([t.loc[m, "imagen"], t.loc[m, "texto"]], [yi, yi],
                color="#BBBBBB", lw=2, zorder=1)
    ax.scatter(t["imagen"], y, s=80, color=COL_IMG, label="Imagen", zorder=3, **EDGE)
    ax.scatter(t["texto"], y, s=80, color=COL_TXT, label="Texto", zorder=3, **EDGE)
    ax.set_yticks(y); ax.set_yticklabels(ORDEN_MODELOS)
    ax.set_xlabel("Tiempo total de ejecución (horas)")
    ax.grid(axis="x", ls=":", alpha=0.5); ax.set_axisbelow(True)
    ax.margins(y=0.2)
    ax.legend(frameon=True, fontsize=9)
    fig.tight_layout()
    guardar(fig, "tiempo_por_modelo")
def main():
    os.makedirs(SALIDA, exist_ok=True)
    r, ram = cargar()
    print("Generando figuras en:", os.path.abspath(SALIDA))

    fig_crossmodal(r)
    for met in ["ARI", "AMI", "NMI"]:
        for mod in ["imagen", "texto"]:
            fig_boxplot(r, met, mod, f"boxplot_{met}_{mod}")
    fig_victorias_externa(r)
    for mod in ["imagen", "texto"]:
        fig_boxplot(r, "Silhouette_mean", mod, f"boxplot_silhouette_{mod}")
    fig_victorias_silhouette(r)
    fig_tiempo(r)
    fig_ram(ram)
    fig_tiempo_por_tamano(r)
    fig_tiempo_por_modelo()

    print("\nListo. Todas las figuras estan en la carpeta '%s' (PNG)." % SALIDA)


if __name__ == "__main__":
    main()