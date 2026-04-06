Este módulo permite descargar datasets desde OpenML y generar subconjuntos de datos para experimentación en algoritmos de clustering.
---
⚙️ Instalación
1. Crear entorno virtual
```bash
python -m venv .venv
```
2. Activar entorno virtual
En Windows (PowerShell):
```bash
.venv\Scripts\Activate.ps1
```
---
3. Instalar dependencias
```bash
pip install -r requirements.txt
```
---
🧩 Configuración

El comportamiento del script se define mediante el archivo:
datasets_config.json
Parámetros

max_rows: número máximo de instancias a usar por dataset

output_dir: carpeta donde se guardarán los datos generados

random_state: semilla para reproducibilidad

groups: conjuntos de datasets organizados por tipo

Ejemplo:
```json
{
  "defaults": {
    "max_rows": 500,
    "output_dir": "data/openml_samples",
    "random_state": 42
  },
  "groups": {
    "1": ["CIFAR-10", "STL-10", "CIFAR-100"],
    "2": ["COIL-100", "USPS", "MNIST", "Fashion-MNIST"],
    "3": ["ImageNet-10", "CUB-200"]
  }
}
```
---
▶️ Uso
Ejecutar por grupo
```bash
python openml_image_sampler.py --group 1
python openml_image_sampler.py --group 2
python openml_image_sampler.py --group 3
```
---
Ejecutar datasets específicos
```bash
python openml_image_sampler.py --datasets CIFAR-10 MNIST Fashion-MNIST
```
---
⚙️ Opciones adicionales
Cambiar número de instancias
```bash
python openml_image_sampler.py --group 1 --max-rows 200
```
Cambiar directorio de salida
```bash
python openml_image_sampler.py --group 3 --output-dir data/pruebas
```
Cambiar semilla
```bash
python openml_image_sampler.py --group 2 --random-state 123
```
Desactivar muestreo estratificado
```bash
python openml_image_sampler.py --group 1 --no-stratify
```
