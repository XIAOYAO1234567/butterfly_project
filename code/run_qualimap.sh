######## BAM Quality Control with Qualimap ########
#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=48:0:0
#$ -l h_vmem=8G
#$ -N run_qualimap_reports

# Starting timestamp
start=$(date +%s)

# Activate qualimap
module load qualimap

# Output directory for cleaned HTML reports
OUT_DIR="mapping_results_qualimap_reports"
mkdir -p "$OUT_DIR"

# Loop over all BAM files
for bam in mapping_results/*.sorted.bam
do
  base=$(basename "$bam" .sorted.bam)
  tmp_dir="${OUT_DIR}/${base}_tmp"

  echo "[INFO] Running Qualimap for $base"

  # Run Qualimap and output to temporary subdirectory
  qualimap bamqc \
    -bam "$bam" \
    -outdir "$tmp_dir" \
    -outformat HTML \
    -nt 4 \
    --java-mem-size=4G

 # Zip the entire HTML report directory
  zip -r "${OUT_DIR}/${base}.qualimap.zip" "$tmp_dir"

  # Clean up temporary directory
  rm -rf "$tmp_dir"

  echo "[INFO] Finished $base"
done

# Finishing timestamp
end=$(date +%s)
runtime=$((end - start))

echo "Finished Qualimap QC"
echo "Runtime: $runtime seconds"
echo $(date)
