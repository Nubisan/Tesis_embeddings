# -*- coding: utf-8 -*-
"""
bench_modelos.py
================
Carga de los 4 modelos (CLIP, BLIP, Qwen2-VL, InternVL) para CPU o GPU.

En CPU: dtype float32 (float16 en CPU o falla o es lentisimo).
En GPU: CLIP/BLIP en float32 y Qwen/InternVL en float16 -> IGUAL que tus
        notebooks de Colab, asi los embeddings coinciden con los canonicos.

Las funciones embed_images/embed_texts aceptan un parametro opcional `batch`
para poder REDUCIRLO en caliente si ocurre CUDA out of memory (el runner lo
usa para reintentar automaticamente con un batch mas pequeno).

cargar_modelo(nombre, device, dtype) -> dict con:
    embed_images(pil_list, batch=None) -> np.ndarray [n, dim]  (float32, sin normalizar)
    embed_texts(text_list, batch=None) -> np.ndarray [n, dim]
    model_name, dtype (str), batch_img, batch_txt
"""

import numpy as np


def _torch_dtype(dtype_str):
    import torch
    return {'float32': torch.float32, 'bfloat16': torch.bfloat16, 'float16': torch.float16}[dtype_str]


# ---------------------------------------------------------------------------
# CLIP
# ---------------------------------------------------------------------------
def _cargar_clip(device, dtype_str, batch_img, batch_txt):
    import torch
    from transformers import CLIPModel, CLIPProcessor
    MODEL_NAME = 'openai/clip-vit-base-patch32'
    model = CLIPModel.from_pretrained(MODEL_NAME).to(device).eval()
    processor = CLIPProcessor.from_pretrained(MODEL_NAME)

    @torch.inference_mode()
    def embed_images(pil_list, batch=None):
        b = batch or batch_img; f = []
        for i in range(0, len(pil_list), b):
            inp = processor(images=pil_list[i:i + b], return_tensors='pt').to(device)
            out = model.get_image_features(pixel_values=inp['pixel_values'])
            emb = out if torch.is_tensor(out) else out.pooler_output
            f.append(emb.float().cpu().numpy())
        return np.concatenate(f, 0)

    @torch.inference_mode()
    def embed_texts(text_list, batch=None):
        b = batch or batch_txt; f = []
        for i in range(0, len(text_list), b):
            inp = processor(text=text_list[i:i + b], return_tensors='pt',
                            padding=True, truncation=True, max_length=77).to(device)
            out = model.get_text_features(input_ids=inp['input_ids'],
                                          attention_mask=inp.get('attention_mask'))
            emb = out if torch.is_tensor(out) else out.pooler_output
            f.append(emb.float().cpu().numpy())
        return np.concatenate(f, 0)

    return MODEL_NAME, embed_images, embed_texts


# ---------------------------------------------------------------------------
# BLIP
# ---------------------------------------------------------------------------
def _cargar_blip(device, dtype_str, batch_img, batch_txt):
    import torch
    from transformers import BlipForImageTextRetrieval, AutoProcessor
    MODEL_NAME = 'Salesforce/blip-itm-base-coco'
    model = BlipForImageTextRetrieval.from_pretrained(MODEL_NAME).to(device).eval()
    processor = AutoProcessor.from_pretrained(MODEL_NAME)

    @torch.inference_mode()
    def embed_images(pil_list, batch=None):
        b = batch or batch_img; f = []
        for i in range(0, len(pil_list), b):
            inp = processor(images=pil_list[i:i + b], return_tensors='pt').to(device)
            vis = model.vision_model(pixel_values=inp['pixel_values'])[0]
            f.append(model.vision_proj(vis[:, 0, :]).float().cpu().numpy())
        return np.concatenate(f, 0)

    @torch.inference_mode()
    def embed_texts(text_list, batch=None):
        b = batch or batch_txt; f = []
        for i in range(0, len(text_list), b):
            inp = processor(text=text_list[i:i + b], return_tensors='pt',
                            padding=True, truncation=True, max_length=512).to(device)
            txt = model.text_encoder(input_ids=inp['input_ids'],
                                     attention_mask=inp['attention_mask'])[0]
            f.append(model.text_proj(txt[:, 0, :]).float().cpu().numpy())
        return np.concatenate(f, 0)

    return MODEL_NAME, embed_images, embed_texts


# ---------------------------------------------------------------------------
# InternVL2-1B
# ---------------------------------------------------------------------------
def _cargar_internvl(device, dtype_str, batch_img, batch_txt):
    import torch
    import torchvision.transforms as T
    from transformers import AutoModel, AutoTokenizer
    MODEL_NAME = 'OpenGVLab/InternVL2-1B'
    dt = _torch_dtype(dtype_str)
    model = AutoModel.from_pretrained(MODEL_NAME, torch_dtype=dt,
                                      trust_remote_code=True).to(device).eval()
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME, trust_remote_code=True, use_fast=False)
    _tf = T.Compose([T.Resize((448, 448)), T.ToTensor(),
                     T.Normalize((0.485, 0.456, 0.406), (0.229, 0.224, 0.225))])

    @torch.inference_mode()
    def embed_images(pil_list, batch=None):
        b = batch or batch_img; f = []
        for i in range(0, len(pil_list), b):
            px = torch.stack([_tf(im.convert('RGB')) for im in pil_list[i:i + b]]).to(device, dtype=dt)
            vit = model.vision_model(pixel_values=px).last_hidden_state
            f.append(vit.mean(1).float().cpu().numpy())
        return np.concatenate(f, 0)

    @torch.inference_mode()
    def embed_texts(text_list, batch=None):
        b = batch or batch_txt; f = []
        for i in range(0, len(text_list), b):
            enc = tokenizer(list(text_list[i:i + b]), return_tensors='pt',
                            padding=True, truncation=True, max_length=512).to(device)
            out = model.language_model(input_ids=enc['input_ids'],
                                       attention_mask=enc['attention_mask'], output_hidden_states=True)
            h = out.hidden_states[-1]
            m = enc['attention_mask'].unsqueeze(-1).to(h.dtype)
            f.append(((h * m).sum(1) / m.sum(1).clamp(min=1)).float().cpu().numpy())
        return np.concatenate(f, 0)

    return MODEL_NAME, embed_images, embed_texts


# ---------------------------------------------------------------------------
# Qwen2-VL-2B
# ---------------------------------------------------------------------------
def _cargar_qwen(device, dtype_str, batch_img, batch_txt):
    import torch
    from transformers import Qwen2VLForConditionalGeneration, AutoProcessor
    from qwen_vl_utils import process_vision_info
    MODEL_NAME = 'Qwen/Qwen2-VL-2B-Instruct'
    dt = _torch_dtype(dtype_str)
    MIN_PIXELS = 4 * 28 * 28
    MAX_PIXELS = 128 * 28 * 28

    model = Qwen2VLForConditionalGeneration.from_pretrained(MODEL_NAME, torch_dtype=dt).to(device).eval()
    model.config.use_cache = False
    processor = AutoProcessor.from_pretrained(MODEL_NAME, min_pixels=MIN_PIXELS, max_pixels=MAX_PIXELS)

    def _mp(h, mask):
        m = mask.unsqueeze(-1).to(h.dtype)
        return (h * m).sum(1) / m.sum(1).clamp(min=1)

    @torch.inference_mode()
    def _embed_batch(inp):
        inp = inp.to(device)
        out = model.model(**inp, use_cache=False, return_dict=True)   # modelo base: sin logits
        emb = _mp(out.last_hidden_state, inp['attention_mask']).float().cpu().numpy()
        del out, inp
        if device == 'cuda':
            torch.cuda.empty_cache()
        return emb

    @torch.inference_mode()
    def embed_texts(text_list, batch=None):
        b = batch or batch_txt; f = []
        for i in range(0, len(text_list), b):
            batch_i = text_list[i:i + b]
            msgs = [[{'role': 'user', 'content': [{'type': 'text', 'text': str(t)}]}] for t in batch_i]
            pr = [processor.apply_chat_template(m, tokenize=False, add_generation_prompt=False) for m in msgs]
            inp = processor(text=pr, padding=True, truncation=True, max_length=512, return_tensors='pt')
            f.append(_embed_batch(inp))
        return np.concatenate(f, 0)

    @torch.inference_mode()
    def embed_images(pil_list, batch=None):
        b = batch or batch_img; f = []
        for i in range(0, len(pil_list), b):
            batch_i = [im.convert('RGB') for im in pil_list[i:i + b]]
            msgs = [[{'role': 'user', 'content': [{'type': 'image', 'image': im},
                                                  {'type': 'text', 'text': '.'}]}] for im in batch_i]
            pr = [processor.apply_chat_template(m, tokenize=False, add_generation_prompt=False) for m in msgs]
            img_in, _ = process_vision_info(msgs)
            inp = processor(text=pr, images=img_in, padding=True, return_tensors='pt')
            f.append(_embed_batch(inp))
        return np.concatenate(f, 0)

    return MODEL_NAME, embed_images, embed_texts


# ---------------------------------------------------------------------------
# Fabrica
# ---------------------------------------------------------------------------
# batch por defecto (los de tus notebooks). En GPU 16 GB puedes subirlos; si un
# dataset de texto largo da CUDA OOM, el runner los baja solo (o usa --batch-txt).
_DEFAULTS = {
    'clip':     dict(cargar=_cargar_clip,     batch_img=64, batch_txt=64, dtype='float32'),
    'blip':     dict(cargar=_cargar_blip,     batch_img=32, batch_txt=32, dtype='float32'),
    'internvl': dict(cargar=_cargar_internvl, batch_img=8,  batch_txt=8,  dtype='float32'),
    'qwen_vl':  dict(cargar=_cargar_qwen,     batch_img=4,  batch_txt=8,  dtype='float32'),
}


def cargar_modelo(nombre, device='cpu', dtype=None, batch_img=None, batch_txt=None):
    if nombre not in _DEFAULTS:
        raise ValueError(f'Modelo desconocido: {nombre}. Opciones: {list(_DEFAULTS)}')
    d = _DEFAULTS[nombre]
    dtype = dtype or d['dtype']
    bi = batch_img or d['batch_img']
    bt = batch_txt or d['batch_txt']
    try:
        model_name, embed_images, embed_texts = d['cargar'](device, dtype, bi, bt)
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            f"Falta una dependencia para '{nombre}': {e}."
        ) from e
    return {
        'modelo': nombre, 'model_name': model_name,
        'embed_images': embed_images, 'embed_texts': embed_texts,
        'dtype': dtype, 'batch_img': bi, 'batch_txt': bt,
    }
