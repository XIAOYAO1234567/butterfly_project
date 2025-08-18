import os, csv, numpy as np, allel

vcf_path     = "bcftools_ml_allchr_10-70/NC_069110.1/merged.NC_069110.1.snps.Q30_DP10-70_MQ30.vcf.gz"
output_dir   = "real_data_npy"
window_size  = 10000
FIXED_ROWS   = 200
FIXED_LENGTH = 600
tsv_out      = os.path.join(output_dir, "windows_index.tsv")

os.makedirs(output_dir, exist_ok=True)

callset   = allel.read_vcf(vcf_path, fields=["variants/POS", "calldata/GT"])
positions = callset["variants/POS"]
genotypes = callset["calldata/GT"]
geno_bin  = (genotypes.sum(axis=2) >= 1).astype(np.int8)

records  = []
basename = os.path.basename(vcf_path).split(".")[1]

start_pos = positions[0] // window_size * window_size
end_pos   = positions[-1]

for win_start in range(start_pos, end_pos, window_size):
    win_end = win_start + window_size
    mask    = (positions >= win_start) & (positions < win_end)
    snp_blk = geno_bin[mask]
    if snp_blk.shape[0] == 0:
        continue
    mat = snp_blk.T
    if mat.shape[0] < FIXED_ROWS:
        reps = (FIXED_ROWS + mat.shape[0] - 1) // mat.shape[0]
        mat  = np.vstack([mat] * reps)[:FIXED_ROWS]
    else:
        mat = mat[:FIXED_ROWS]
    if mat.shape[1] < FIXED_LENGTH:
        mat = np.pad(mat, ((0, 0), (0, FIXED_LENGTH - mat.shape[1])), mode="constant")
    else:
        mat = mat[:, :FIXED_LENGTH]
    out_name = f"{basename}_{win_start}_{win_end}.npy"
    out_path = os.path.join(output_dir, out_name)
    np.save(out_path, mat)
    records.append([basename, win_start, win_end, out_path, snp_blk.shape[0]])

with open(tsv_out, "w", newline="") as f:
    csv.writer(f, delimiter="\t").writerows(
        [["chromosome", "start", "end", "npy_path", "num_snps"]] + records
    )

