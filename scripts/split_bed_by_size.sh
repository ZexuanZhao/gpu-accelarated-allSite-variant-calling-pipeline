#!/bin/bash

input_bed="$1"
output_prefix="$2"
max_total_size="$3"  # Maximum total region size

current_total_size=0
file_index=1
exec 3> "${output_prefix}.${file_index}.bed"

while IFS=$'\t' read -r chrom start end rest; do
  region_size=$((end - start))

  if [[ $((current_total_size + region_size)) -le "$max_total_size" ]]; then
    echo -e "$chrom\t$start\t$end\t$rest" >&3
    current_total_size=$((current_total_size + region_size))
  else
    exec 3>&-
    file_index=$((file_index + 1))
    exec 3> "${output_prefix}.${file_index}.bed"
    echo -e "$chrom\t$start\t$end\t$rest" >&3
    current_total_size=$region_size
  fi
done < "$input_bed"

exec 3>&-
#echo "Split BED file into files with maximum total region size of $max_total_size."