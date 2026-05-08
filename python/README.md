# Clustering Python

Migración a Python de los algoritmos de clustering originalmente en R.

## Instalación

Con Python 3.11 ya instalado, desde la carpeta del proyecto:

```bash
py -3.11 -m venv venv
venv\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt
```

## Estructura

```
clustering-python/
├── algorithms/             # módulos de cada algoritmo (1 por archivo)
│   └── __init__.py
├── datasets_local/         # caché de datasets descargados de OpenML
├── openml_loader.py        # equivalente a Openml.R
├── testing.py              # equivalente a Testing.R (ejecutor)
├── peakRAM_log.csv         # log de tiempos y memoria pico (se crea al ejecutar)
├── requirements.txt
└── README.md
```

## Uso

### Probar el loader (descarga datasets):
```bash
python openml_loader.py
```

### Ejecutar un algoritmo individual:
```bash
python testing.py 1     # 1 = Bat, 2 = CSCLP, ... 8 = SCK1
```

### Ejecutar todos en paralelo:
```bash
python testing.py 9
```

## Mapeo de algoritmos

| numalt | Algoritmo | Módulo Python                |
|--------|-----------|------------------------------|
| 1      | Bat       | `algorithms/bat.py`          |
| 2      | CSCLP     | `algorithms/csclp.py`        |
| 3      | Kmedoids  | `algorithms/kmedoids.py`     |
| 4      | ACO       | `algorithms/aco.py`          |
| 5      | PSO       | `algorithms/pso.py`          |
| 6      | HCAKC     | `algorithms/hcakc.py`        |
| 7      | KM-MILP   | `algorithms/km_milp.py`      |
| 8      | SCK1      | `algorithms/sck1_final.py`   |

Cada módulo debe exponer una función:

```python
def run(odatasets_unique):
    # odatasets_unique es un pandas.DataFrame con columna 'dataset' que
    # contiene cada dataset como DataFrame
    ...
```
