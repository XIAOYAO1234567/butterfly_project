#!/bin/bash

REGION_FILE="ld_candidate_regions.tsv"
INPUT_DIR="ngsld_output"
OUTPUT_DIR="ld_tables_candidate"
mkdir -p "$OUTPUT_DIR"

while IFS=$'\t' read -r chr start end || [[ -n "$chr" ]]; do
    chr=$(echo "$chr" | xargs)
    start=$(echo "$start" | xargs)
    end=$(echo "$end" | xargs)

    [[ -z "$chr" || -z "$start" || -z "$end" ]] && continue

    found=0
    chr_regex=$(echo "$chr" | sed 's/\./\\./g')  

    for f in ${INPUT_DIR}/ngsld_${chr}_*.filtered.txt.gz; do
        fname=$(basename "$f")

        if [[ "$fname" =~ ngsld_${chr_regex}_([0-9]+)_([0-9]+)\.filtered\.txt\.gz ]]; then
            window_start=${BASH_REMATCH[1]}
            window_end=${BASH_REMATCH[2]}

            if (( start <= window_end && end >= window_start )); then
                base=$(basename "$f" .filtered.txt.gz)
                outfile="${OUTPUT_DIR}/${base}_${start}_${end}.ld.tsv"

                echo " Extracting $chr:$start-$end from $fname"
                zcat "$f" | awk -v s=$start -v e=$end -v chr=$chr '
                    BEGIN { FS="\t"; OFS="\t" }
                    NR > 1 {
                        split($1, a, ":");
                        split($2, b, ":");
                        pos1 = a[2];
                        pos2 = b[2];
                        if (a[1]==chr && b[1]==chr && pos1 >= s && pos1 <= e && pos2 >= s && pos2 <= e) {
                            print pos1, pos2, $7;
                        }
                    }
                ' > "$outfile"
                found=1
                break
            fi
        fi
    done

    if [[ $found -eq 0 ]]; then
        echo " No matching file found for $chr $start $end"
    fi
done < "$REGION_FILE"
