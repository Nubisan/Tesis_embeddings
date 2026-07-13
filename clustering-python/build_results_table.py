# -*- coding: utf-8 -*-
"""
Combina los archivos pred_*.csv de cada modelo y modalidad en un unico CSV
en formato largo (una fila por dataset x algoritmo x modelo x modalidad).

Ejecutar:  python build_results_table.py
Genera 'resultados_largo.csv' en la carpeta de salida.
"""

import os
import glob
import numpy as np
import pandas as pd

# ====================
# 1) CONFIGURACION
# ===================

# Mapeo (modelo, modalidad) -> carpeta que contiene los pred_*.csv.
FUENTES = {
    ("CLIP",      "imagen"): r"C:\Users\DELL\Documents\8semestre\Tesis\Tesis_embeddings\clustering-python\predictions\clip\imagen",
    ("BLIP",      "imagen"): r"C:\Users\DELL\Documents\8semestre\Tesis\Tesis_embeddings\clustering-python\predictions\blip\imagen",
    ("Qwen2-VL",  "imagen"): r"C:\Users\DELL\Documents\8semestre\Tesis\Tesis_embeddings\clustering-python\predictions\qwen_vl\imagen",
    ("InternVL2", "imagen"): r"C:\Users\DELL\Documents\8semestre\Tesis\Tesis_embeddings\clustering-python\predictions\internvl\imagen",
    ("CLIP",      "texto"):  r"C:\Users\DELL\Documents\8semestre\Tesis\Tesis_embeddings\clustering-python\predictions\clip\texto",
    ("BLIP",      "texto"):  r"C:\Users\DELL\Documents\8semestre\Tesis\Tesis_embeddings\clustering-python\predictions\blip\texto",
    ("Qwen2-VL",  "texto"):  r"C:\Users\DELL\Documents\8semestre\Tesis\Tesis_embeddings\clustering-python\predictions\qwen_vl\texto",
    ("InternVL2", "texto"):  r"C:\Users\DELL\Documents\8semestre\Tesis\Tesis_embeddings\clustering-python\predictions\internvl\texto",
}

SALIDA = "resultados_largo.csv"

# ============================================================
# 2) Mapeo nombre de archivo -> nombre del algoritmo
# ============================================================
ALGO_MAP = {
    "pred_ACO.csv":      "K-MeansACO",
    "pred_BAT.csv":      "K-MeansBA",
    "pred_CSCLP.csv":    "CSCLP",
    "pred_HCAKC.csv":    "HCAKC",
    "pred_KMEDOIDS.csv": "K-MedoidsSC",
    "pred_KM-MILP.csv":  "MILP-KM",
    "pred_PSO.csv":      "K-MeansPSO",
    "pred_SCK1.csv":     "SCK",
}

# Columnas de metricas que queremos conservar.
# Se descartan y_predict / y_true.
COLS_METRICAS = [
    "name", "n", "k",
    "target_cardinality", "cardinality_pred",
    "Execution_Time", "ARI", "AMI", "NMI", "Silhouette_mean",
]


def error_cardinalidad(row):
    """Error % entre cardinalidad objetivo y predicha, comparando como
    CONJUNTO ordenado """
    try:
        t = np.sort(np.array(str(row["target_cardinality"]).split(), dtype=float))
        p = np.sort(np.array(str(row["cardinality_pred"]).split(), dtype=float))
    except Exception:
        return np.nan
    if len(t) != len(p) or t.sum() == 0:
        return np.nan
    return np.abs(t - p).sum() / t.sum() * 100.0


def cargar_carpeta(carpeta, modelo, modalidad):
    """Lee los pred_*.csv de una carpeta y devuelve un DataFrame largo."""
    filas = []
    for archivo, algo in ALGO_MAP.items():
        ruta = os.path.join(carpeta, archivo)
        if not os.path.exists(ruta):
            print(f"   [aviso] falta {archivo} en {carpeta}")
            continue
        # Leer solo las columnas de metricas que existan en el csv
        cabecera = pd.read_csv(ruta, nrows=0).columns.tolist()
        usar = [c for c in COLS_METRICAS if c in cabecera]
        d = pd.read_csv(ruta, usecols=usar)
        d["algoritmo"] = algo
        d["modelo"] = modelo
        d["modalidad"] = modalidad
        filas.append(d)
    if not filas:
        return None
    return pd.concat(filas, ignore_index=True)


def main():
    bloques = []
    for (modelo, modalidad), carpeta in FUENTES.items():
        if not os.path.isdir(carpeta):
            print(f"[saltado] {modelo}/{modalidad}: no existe la carpeta '{carpeta}'")
            continue
        print(f"[leyendo] {modelo}/{modalidad}  ({carpeta})")
        df = cargar_carpeta(carpeta, modelo, modalidad)
        if df is not None:
            bloques.append(df)

    if not bloques:
        print("\nNo se leyo ningun archivo. Revisa las rutas en FUENTES.")
        return

    data = pd.concat(bloques, ignore_index=True)

    # Error de cardinalidad (conjunto ordenado)
    if {"target_cardinality", "cardinality_pred"}.issubset(data.columns):
        data["Error_Card_Pct"] = data.apply(error_cardinalidad, axis=1)

    # Orden de columnas amigable
    primero = ["name", "algoritmo", "modelo", "modalidad", "n", "k"]
    resto = [c for c in data.columns if c not in primero]
    data = data[[c for c in primero if c in data.columns] + resto]

    data.to_csv(SALIDA, index=False, encoding="utf-8-sig")

    print("\n================= RESUMEN =================")
    print(f"Archivo generado : {os.path.abspath(SALIDA)}")
    print(f"Filas totales    : {len(data)}")
    print(f"Modelos          : {sorted(data['modelo'].unique())}")
    print(f"Modalidades      : {sorted(data['modalidad'].unique())}")
    print(f"Algoritmos       : {data['algoritmo'].nunique()}")
    print(f"Datasets         : {data['name'].nunique()}")
    print("\nFilas por modelo x modalidad:")
    print(data.groupby(['modelo', 'modalidad']).size().to_string())


if __name__ == "__main__":
    main()