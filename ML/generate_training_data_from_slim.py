#!/usr/bin/env python3
import os, re, numpy as np
from glob import glob
from pathlib import Path

ROOT_DIR = "SLiM_update/outputs"
OUTPUT_DIR = "flexsweep_training_data_200hap_100kb"
os.makedirs(OUTPUT_DIR, exist_ok=True)

WINDOW_SIZE = 10_000
SAMPLE_ROWS = 200
np.random.seed(42)

genome_lengths = {
    "NC_069085.1": 19462324,
    "NC_069087.1": 19123075,
    "NC_069090.1": 17880892,
    "NC_069093.1": 17486486,
    "NC_069101.1": 15486372,
    "NC_069105.1": 13224086,
    "NC_069110.1": 21498244
}

selection_regions = {
    "NC_069110.1": {"balancing": [(17400001, 17420000), (5920001, 5945000)]},
    "NC_069105.1": {"sweep": [(5050001, 5070000)], "neutral": [(5150001, 5200000)]},
    "NC_069101.1": {"sweep": [(6810001, 6840000)], "neutral": [(6740001, 6810000)]},
    "NC_069093.1": {"neutral": [(1070001, 1210000)], "balancing": [(1210001, 1230000)]},
    "NC_069090.1": {"sweep": [(11430001, 11470000)], "balancing": [(1360001, 1380000)], "neutral": [(1380001, 1420000), (11403001, 11430000)]},
    "NC_069087.1": {"balancing": [(890001, 910000), (1110001, 1130000), (1650001, 1680000)], "neutral": [(4100001, 4400000)]},
    "NC_069085.1": {"sweep": [(10725001, 10750000), (11040001, 11060000)], "neutral": [(10860001, 10920000)]}
}

def assign_label(chrom, mid):
    for lbl in ("balancing", "sweep", "neutral"):
        for s, e in selection_regions.get(chrom, {}).get(lbl, []):
            if s <= mid <= e:
                return lbl
    return None

def parse_ms(path):
    with open(path) as f:
        lines = f.readlines()
    pos_line = next(l for l in lines if l.startswith("positions:"))
    pos = np.array(pos_line.split()[1:], dtype=float)
    seq = [l.strip() for l in lines if re.fullmatch(r"[01]+", l.strip())]
    gen = np.array([[int(ch) for ch in s] for s in seq], dtype=np.int8)
    return pos, gen

def main():
    recs = []
    for chrom, g_len in genome_lengths.items():
        chrom_dir = os.path.join(ROOT_DIR, chrom)
        pattern = os.path.join(chrom_dir, "replicate", "rep*", f"{chrom}_1000_generations.ms")
        for ms_path in sorted(glob(pattern)):
            pos, gen = parse_ms(ms_path)
            if pos.size == 0 or gen.shape[0] < SAMPLE_ROWS:
                continue
            abs_pos = (pos * g_len).astype(int)
            for start in range(1, g_len, WINDOW_SIZE):
                end = start + WINDOW_SIZE - 1
                mid = (start + end) // 2
                idx = np.where((abs_pos >= start) & (abs_pos <= end))[0]
                if idx.size < 2:
                    continue
                sub = gen[:, idx]
                if sub.shape[0] < SAMPLE_ROWS:
                    continue
                sub = sub[np.random.choice(sub.shape[0], SAMPLE_ROWS, replace=False)]
                lbl = assign_label(chrom, mid)
                if lbl is None:
                    continue
                rep = Path(ms_path).parent.name
                out_path = os.path.join(OUTPUT_DIR, f"{chrom}_{rep}_{start}_{end}.npy")
                np.save(out_path, sub)
                recs.append((out_path, lbl))
    with open(os.path.join(OUTPUT_DIR, "train.tsv"), "w") as fw:
        fw.write("filepath\tlabel\n")
        for fp, lb in recs:
            fw.write(f"{fp}\t{lb}\n")
    print(f"Saved {len(recs)} samples to {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
