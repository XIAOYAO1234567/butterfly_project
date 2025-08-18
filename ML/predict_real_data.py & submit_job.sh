#!/usr/bin/env python3
import torch
import numpy as np
import pandas as pd
from cnn_model import CNN5

npy_index_file = "real_data_npy/windows_index.tsv"
output_tsv     = "real_data_npy/prediction_results.tsv"
model_path     = "best_model_200row.pt"
L              = 600
KEEP_ROWS      = 200
device         = torch.device("cuda" if torch.cuda.is_available() else "cpu")
label_map      = {0: "balancing", 1: "neutral", 2: "sweep"}

df = pd.read_csv(npy_index_file, sep="\t")

model = CNN5(L=L, n_cls=3, keep_rows=KEEP_ROWS).to(device)
model.load_state_dict(torch.load(model_path, map_location=device))
model.eval()

@torch.no_grad()
def predict(path):
    x = np.load(path)
    x = np.pad(x, ((0, 0), (0, max(0, L - x.shape[1]))), mode="constant")[:, :L]
    r = x.shape[0]
    if r < KEEP_ROWS:
        x = np.vstack([x, np.zeros((KEEP_ROWS - r, L), dtype=np.int8)])
    elif r > KEEP_ROWS:
        x = x[:KEEP_ROWS]
    t = torch.tensor(x, dtype=torch.float32).unsqueeze(0).to(device)
    p = torch.softmax(model(t), dim=1).cpu().numpy().squeeze()
    return int(p.argmax()), p

preds, probs = [], []
for p in df["npy_path"]:
    y, pr = predict(p)
    preds.append(label_map[y])
    probs.append(pr)

df["predicted_label"] = preds
df[[f"prob_{label_map[i]}" for i in range(3)]] = pd.DataFrame(probs)
df.to_csv(output_tsv, sep="\t", index=False)
print(f"results written to {output_tsv}")



#!/bin/bash
#$ -cwd
#$ -l h_vmem=10G
#$ -pe smp 4

source ~/miniforge3/etc/profile.d/conda.sh
conda activate flexsweep_env

python predict_real_data.py
