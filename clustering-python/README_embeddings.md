# Embeddings Pipeline - Guía Rápida

Esta guía describe el flujo para correr los 8 algoritmos de clustering sobre embeddings generados con modelos pre-entrenados (CLIP, BLIP, Qwen2-VL, InternVL2).

## Flujo general

```
[1] Colab (GPU)              [2] Tu PC local
─────────────────             ─────────────────
Genera embeddings    ─.zip─>  Descomprime en embeddings/
con CLIP/BLIP/etc.            Corre python testing.py 9 --model clip
                              Revisa CSVs en predictions/
```

---

## Estructura del proyecto

```
clustering-python/
├── colab/
│   └── 01_clip_embeddings.ipynb       # Notebook para Colab
├── embeddings/                         # Vacío al inicio
│   └── clip/                           # Se rellena después de Colab
│       ├── bbc_sport.npz
│       ├── bbc_news.npz
│       └── ...
├── algorithms/                         # Los 8 algoritmos (sin cambios)
├── predictions/                        # Output de cada algoritmo
├── embeddings_loader.py                # Carga .npz (nuevo)
├── testing.py                          # Orquestador (actualizado)
└── README_embeddings.md                # Esta guía
```

---

## Paso 1: Generar embeddings en Colab

1. Subir el notebook `colab/01_clip_embeddings.ipynb` a Google Colab
2. **Runtime → Change runtime type → GPU (T4)**
3. **Runtime → Run all** (ejecuta todo de corrido)
4. Cuando termine (~5-15 min), descargar `clip_embeddings.zip`

**8 datasets generados:**

| Dataset | Tipo | Clases | N |
|---|---|---|---|
| bbc_sport | Texto | 5 | 737 |
| bbc_news | Texto | 5 | 2,225 |
| 20ng_2classes | Texto | 2 | ~1,000 |
| 20ng_3classes | Texto | 3 | ~1,500 |
| ag_news_4k | Texto | 4 | 4,000 |
| orl_faces | Imagen | 40 | 400 |
| mnist_5k | Imagen | 10 | 5,000 |
| fashion_mnist_5k | Imagen | 10 | 5,000 |

---

## Paso 2: Descomprimir en tu PC

1. Mover `clip_embeddings.zip` a la raíz del proyecto
2. Descomprimir:
   ```cmd
   tar -xf clip_embeddings.zip -C embeddings\
   ```
   o con cualquier programa de unzip (Windows lo hace nativo con clic derecho → Extraer)

3. Verificar la estructura:
   ```
   embeddings/clip/
   ├── bbc_sport.npz
   ├── bbc_news.npz
   ├── 20ng_2classes.npz
   ├── 20ng_3classes.npz
   ├── ag_news_4k.npz
   ├── orl_faces.npz
   ├── mnist_5k.npz
   └── fashion_mnist_5k.npz
   ```

4. (Opcional) Inspeccionar contenido:
   ```cmd
   venv\Scripts\activate
   python embeddings_loader.py clip
   ```

   Debería imprimir un resumen con los 8 datasets.

---

## Paso 3: Ejecutar algoritmos

### Modo individual (un algoritmo a la vez)

```cmd
venv\Scripts\activate

python testing.py 1 --model clip    # BAT
python testing.py 2 --model clip    # CSCLP
python testing.py 3 --model clip    # K-Medoids
python testing.py 4 --model clip    # ACO
python testing.py 5 --model clip    # PSO
python testing.py 6 --model clip    # HCAKC
python testing.py 7 --model clip    # KM-MILP
python testing.py 8 --model clip    # SCK1
```

### Modo paralelo (los 8 a la vez, recomendado para terminar rápido)

```cmd
python testing.py 9 --model clip
```

Esto lanza 8 procesos paralelos (uno por algoritmo). Tarda lo que tarde el más lento (típicamente ACO o BAT).

### Resultados

Cada algoritmo escribe un CSV en `predictions/`:
```
predictions/
├── pred_BAT.csv
├── pred_CSCLP.csv
├── pred_KMEDOIDS.csv
├── pred_ACO.csv
├── pred_PSO.csv
├── pred_HCAKC.csv
├── pred_KM-MILP.csv
└── pred_SCK1.csv
```

Cada CSV tiene una fila por dataset con: `name`, `n`, `k`, `target_cardinality`, `cardinality_pred`, `Execution_Time`, `ARI`, `AMI`, `NMI`, `Silhouette_mean`.

---

## Paso 4 (próximo): Otros modelos

Cuando quieras añadir BLIP, Qwen2-VL o InternVL2:

1. Generar notebooks `02_blip_embeddings.ipynb`, etc. (mismo patrón que el de CLIP)
2. Ejecutar en Colab
3. Descomprimir en `embeddings/blip/`, `embeddings/qwen_vl/`, `embeddings/internvl/`
4. Correr con `--model blip`, `--model qwen_vl`, `--model internvl`

Los 8 algoritmos **no requieren cambios** para soportar nuevos modelos. Lo único que cambia es la fuente de embeddings.

---

## Troubleshooting

### Error: "No se encontraron archivos .npz"
- Verifica que existe la carpeta `embeddings/clip/` con los `.npz` dentro
- No deben estar en una subcarpeta extra (ej: `embeddings/clip/clip/`)

### Error: "Out of memory" en algún algoritmo
- Datasets >10k pueden requerir mucha RAM por la matriz N×N
- Empieza con datasets más pequeños y ve escalando
- Si es persistente, contactar para añadir modo chunked

### Algoritmo tarda mucho (>10 min en datasets pequeños)
- ACO es naturalmente lento (50 hormigas × 20 iter × kmeans interno)
- CSCLP/KM-MILP/SCK1 dependen del solver MILP de scipy
- BAT, PSO, HCAKC, K-Medoids son rápidos

---

## Notas técnicas

- **Normalización L2**: los embeddings se guardan con `||x||₂ = 1`, lo que hace que `dist_coseno(u,v) = 1 − u·v`. Esto es mucho más rápido que `pdist(metric="cosine")` clásico.
- **Formato float32**: ocupa la mitad que float64 sin pérdida práctica de precisión.
- **Compatibilidad**: los 8 algoritmos no detectan que los datos son embeddings. Funcionan igual que con datos OpenML tabulares, gracias al "shortcut embeddings" en `_common.py::prepare_data`.
