#!/usr/bin/env python3
import numpy as np, pandas as pd, torch, matplotlib.pyplot as plt
import torch.nn.functional as F
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, ConfusionMatrixDisplay, roc_curve, auc
from sklearn.preprocessing import LabelEncoder, label_binarize
from torch.utils.data import Dataset, DataLoader, Subset
from cnn_model import CNN5

INPUT_TSV     = "flexsweep_training_data_200hap_100kb/train.tsv"
BEST_PTH      = "best_model_200row.pt"
KEEP_ROWS     = 200
FIXED_LENGTH  = 600
VAL_RATIO     = 0.2
BATCH_SIZE    = 32
SEED          = 42
NUM_CLASSES   = 3
DEVICE        = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class SNPWindowDataset(Dataset):
    def __init__(self, tsv, keep_rows, fixed_len):
        self.paths, self.labels, self.keep, self.L = [], [], keep_rows, fixed_len
        with open(tsv) as f:
            next(f)
            for line in f:
                p, lb = line.strip().split("\t")
                self.paths.append(p)
                self.labels.append(lb)
        self.le = LabelEncoder()
        self.y  = self.le.fit_transform(self.labels)
        self.rng = np.random.RandomState(SEED)

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):
        x = np.load(self.paths[idx])
        r, c = x.shape
        if r >= self.keep:
            x = x[self.rng.choice(r, self.keep, replace=False)]
        else:
            pad = self.keep - r
            x = np.vstack([x, np.zeros((pad, c), dtype=x.dtype)])
        if c > self.L:
            x = x[:, : self.L]
        elif c < self.L:
            x = np.pad(x, ((0, 0), (0, self.L - c)), mode="constant")
        return torch.tensor(x, dtype=torch.float32), torch.tensor(self.y[idx], dtype=torch.long)

ds = SNPWindowDataset(INPUT_TSV, KEEP_ROWS, FIXED_LENGTH)
idx_tr, idx_te = train_test_split(np.arange(len(ds)), test_size=VAL_RATIO, random_state=SEED, stratify=ds.y)
loader = DataLoader(Subset(ds, idx_te), batch_size=BATCH_SIZE)

model = CNN5(L=FIXED_LENGTH, n_cls=NUM_CLASSES, keep_rows=KEEP_ROWS).to(DEVICE)
model.load_state_dict(torch.load(BEST_PTH, map_location=DEVICE))
model.eval()

y_true, y_pred, y_prob = [], [], []
with torch.no_grad():
    for x, y in loader:
        x = x.to(DEVICE)
        prob = F.softmax(model(x), dim=1).cpu().numpy()
        y_true.extend(y.numpy())
        y_pred.extend(prob.argmax(axis=1))
        y_prob.extend(prob)

y_true = np.array(y_true)
y_pred = np.array(y_pred)
y_prob = np.array(y_prob)
labels = ds.le.classes_

rep = classification_report(y_true, y_pred, target_names=labels, output_dict=True)
pd.DataFrame(rep).T.to_csv("classification_report.csv")

cm = confusion_matrix(y_true, y_pred)
ConfusionMatrixDisplay(cm, display_labels=labels).plot(cmap="Blues", values_format="d")
plt.title("Confusion Matrix")
plt.savefig("confusion_matrix_test.png")
plt.close()

y_bin = label_binarize(y_true, classes=np.arange(NUM_CLASSES))
fpr, tpr, roc_auc = {}, {}, {}
for i in range(NUM_CLASSES):
    fpr[i], tpr[i], _ = roc_curve(y_bin[:, i], y_prob[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])
plt.figure(figsize=(6, 5))
for i in range(NUM_CLASSES):
    plt.plot(fpr[i], tpr[i], label=f"{labels[i]} AUC={roc_auc[i]:.2f}")
plt.plot([0, 1], [0, 1], "k--")
plt.xlabel("FPR")
plt.ylabel("TPR")
plt.title("ROC")
plt.legend()
plt.tight_layout()
plt.savefig("roc_curve_test.png")
plt.close()

pd.DataFrame({
    "true": [labels[i] for i in y_true],
    "pred": [labels[i] for i in y_pred]
}).to_csv("prediction_results.csv", index=False)

print("evaluation finished")
