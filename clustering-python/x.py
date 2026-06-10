"""
Script para explorar la estructura del archivo NPZ de los embeddings de CLIP.
"""

import numpy as np

npz = np.load('embeddings/internvl/20ng_2classes.npz', allow_pickle=True)
print(list(npz.files))        # Claves
print(npz['X'].shape)         # Ej: (400, 512)
print(npz['y'][:10])          # Primeros 10 labels
print(npz['class_names'])     # Nombres de clases
print(npz['metadata'][0])     # Metadata
