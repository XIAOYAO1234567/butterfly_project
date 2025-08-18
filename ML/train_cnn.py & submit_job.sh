#!/usr/bin/env python3
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, torch, torch.nn as nn, torch.optim as optim, torch.nn.functional as F
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import Dataset, DataLoader, Subset

INPUT_TSV    = "flexsweep_training_data_200hap_100kb/train.tsv"
KEEP_ROWS    = 200
FIXED_LENGTH = 600
VAL_RATIO    = 0.2
BATCH_SIZE   = 32
EPOCHS       = 50
LR           = 1e-3
PATIENCE     = 5
SEED         = 42
NUM_CLASSES  = 3
device       = torch.device("cuda" if torch.cuda.is_available() else "cpu")
rng          = np.random.RandomState(SEED)

class SNPWindowDataset(Dataset):
    def __init__(self, tsv, keep, length):
        self.paths, self.labels = [], []
        with open(tsv) as f:
            next(f)
            for line in f:
                p, lb = line.strip().split("\t")
                self.paths.append(p)
                self.labels.append(lb)
        self.keep, self.L = keep, length
        self.le = LabelEncoder()
        self.y = self.le.fit_transform(self.labels)
    def __len__(self):
        return len(self.paths)
    def __getitem__(self, idx):
        x = np.load(self.paths[idx])
        r, c = x.shape
        if r >= self.keep:
            x = x[rng.choice(r, self.keep, replace=False)]
        else:
            x = np.vstack([x, np.zeros((self.keep - r, c), np.int8)])
        x = x[:, :self.L] if c >= self.L else np.pad(x, ((0, 0), (0, self.L - c)))
        return torch.tensor(x, dtype=torch.float32), torch.tensor(self.y[idx])

class CNN5(nn.Module):
    def __init__(self, L=FIXED_LENGTH, n_cls=NUM_CLASSES, keep_rows=KEEP_ROWS):
        super().__init__()
        self.conv1 = nn.Conv1d(keep_rows, 64, 5, padding=2)
        self.pool1 = nn.MaxPool1d(2)
        self.conv2 = nn.Conv1d(64, 128, 3, padding=1)
        self.pool2 = nn.MaxPool1d(2)
        self.dp    = nn.Dropout(0.5)
        self.fc1   = nn.Linear(128 * (L // 4), 128)
        self.fc2   = nn.Linear(128, n_cls)
    def forward(self, x):
        x = self.pool1(F.relu(self.conv1(x)))
        x = self.pool2(F.relu(self.conv2(x)))
        x = x.view(x.size(0), -1)
        x = self.dp(x)
        x = F.relu(self.fc1(x))
        return self.fc2(x)

ds = SNPWindowDataset(INPUT_TSV, KEEP_ROWS, FIXED_LENGTH)
tr_idx, va_idx = train_test_split(np.arange(len(ds)), test_size=VAL_RATIO, random_state=SEED, stratify=ds.y)
tr_loader = DataLoader(Subset(ds, tr_idx), batch_size=BATCH_SIZE, shuffle=True)
va_loader = DataLoader(Subset(ds, va_idx), batch_size=BATCH_SIZE)

criterion = nn.CrossEntropyLoss().to(device)
model = CNN5().to(device)
opt = optim.Adam(model.parameters(), lr=LR)

tr_loss_hist, va_loss_hist, tr_acc_hist, va_acc_hist = [], [], [], []
best_acc, wait = 0, 0

for epoch in range(1, EPOCHS + 1):
    model.train(); tr_corr = tr_loss = 0
    for x, y in tr_loader:
        x, y = x.to(device), y.to(device)
        opt.zero_grad()
        out = model(x)
        loss = criterion(out, y)
        loss.backward(); opt.step()
        tr_loss += loss.item() * x.size(0)
        tr_corr += (out.argmax(1) == y).sum().item()
    tr_loss /= len(tr_idx)
    tr_acc = tr_corr / len(tr_idx)
    model.eval(); va_corr = va_loss = 0
    with torch.no_grad():
        for x, y in va_loader:
            x, y = x.to(device), y.to(device)
            out = model(x)
            va_loss += criterion(out, y).item() * x.size(0)
            va_corr += (out.argmax(1) == y).sum().item()
    va_loss /= len(va_idx)
    va_acc = va_corr / len(va_idx)
    tr_loss_hist.append(tr_loss); va_loss_hist.append(va_loss)
    tr_acc_hist.append(tr_acc);   va_acc_hist.append(va_acc)
    print(f"epoch {epoch:02d}  train_acc {tr_acc:.4f}  val_acc {va_acc:.4f}")
    if va_acc > best_acc:
        best_acc, wait = va_acc, 0
        torch.save(model.state_dict(), "best_model_200row.pt")
    else:
        wait += 1
        if wait >= PATIENCE:
            break

epochs = np.arange(1, len(tr_loss_hist) + 1)
plt.figure(figsize=(8, 4))
plt.subplot(1, 2, 1)
plt.plot(epochs, tr_loss_hist, label="train")
plt.plot(epochs, va_loss_hist, label="val")
plt.xlabel("epoch"); plt.ylabel("loss"); plt.legend()
plt.subplot(1, 2, 2)
plt.plot(epochs, tr_acc_hist, label="train")
plt.plot(epochs, va_acc_hist, label="val")
plt.xlabel("epoch"); plt.ylabel("accuracy"); plt.legend()
plt.tight_layout()
plt.savefig("learning_curves.png")

      
#!/bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l gpu=1               
#$ -pe smp 8             

source ~/.bashrc
conda activate flexsweep_env

python train_cnn_model.py
