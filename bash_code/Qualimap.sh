#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=48:0:0
#$ -l h_vmem=8G

start=$(date +%s)

module load qualimap
OUT_DIR="mapping_results_qualimap_reports"
mkdir -p "$OUT_DIR"

for bam in mapping_results/*.sorted.bam; do
  base=$(basename "$bam" .sorted.bam)
  tmp_dir="${OUT_DIR}/${base}_tmp"
  qualimap bamqc -bam "$bam" -outdir "$tmp_dir" -outformat HTML -nt 4 --java-mem-size=4G
  zip -r "${OUT_DIR}/${base}.qualimap.zip" "$tmp_dir"
  rm -rf "$tmp_dir"
done

end=$(date +%s)
runtime=$((end - start))

echo "Finished Qualimap QC"
echo "Runtime: $runtime seconds"
echo "$(date)"
