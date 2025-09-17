#!/bin/bash
set -euo pipefail

# GATA1 (chrX)
# MYC (chr8q24)
# RAE1 (chr20) RNA export factor

# --------------------------
# Config
# --------------------------
input_bed_filename="hg38_Masked.bed"                      		# the tiny downloadable BED
output_bed_filename="hg38_Masked_Sorted.bed"              		# the one you'll load in IGV
input_bam_filename="K562_erythroleukemia_H3K4me3_sorted.bam"  	# must be coordinate-sorted
demo_region="chr20:57349262-57355972"                     		# your masked locus (also appended below)
demo_bam="demo.bam"
demo_sam="demo.sam"
spike_sam="spike.sam"                                     		# output of your Python spiker
spike_bam="spike.bam"
spike_bam_sorted="spike.sorted.bam"
merged_bam="K562_erythroleukemia_H3K4me3_merged.bam"
merged_bam_sorted="${input_bam_filename}"                 		# intentionally overwrite input_bam

echo "Preparing masked BED (keeping sorted BED for IGV)..."
# Impute/append the demo region (idempotent append may duplicate if re-run; optional 'grep -qxF' guard)
echo -e "chr20\t57349262\t57355972" >> "${input_bed_filename}"
sort -k1,1 -k2,2n "${input_bed_filename}" | uniq > "${output_bed_filename}"

echo "==> Slicing a tiny demo window from the big BAM..."
samtools view -b "${input_bam_filename}" "${demo_region}" -o "${demo_bam}"
samtools index "${demo_bam}"
samtools idxstats "${demo_bam}" | head -n 5
samtools flagstat "${demo_bam}"

echo "==> Converting slice to SAM for spiking..."
samtools view -h "${demo_bam}" > "${demo_sam}"

echo "==> Fabricating duplicate pile (Python)..."
# Your Python script should read demo.sam and write spike.sam
python _script.py   # make sure it writes '${spike_sam}'

echo "==> Converting spiked SAM -> BAM and sorting..."
samtools view -bS "${spike_sam}" | samtools sort -o "${spike_bam_sorted}"
samtools index "${spike_bam_sorted}"

echo "==> (Optional) sort the slice as well to be safe for merge..."
samtools sort -o "${demo_bam%.bam}.sorted.bam" "${demo_bam}"
mv "${demo_bam%.bam}.sorted.bam" "${demo_bam}"
samtools index "${demo_bam}"

echo "==> Merging original BAM + spiked demo BAM..."
# Use -f to overwrite if exists; both inputs coordinate-sorted
samtools merge -f -@8 "${merged_bam}" "${input_bam_filename}" "${spike_bam_sorted}"

echo "==> Sorting and indexing the merged BAM..."
samtools sort -@8 -o "${merged_bam_sorted}" "${merged_bam}"
samtools index "${merged_bam_sorted}"

echo "==> Cleaning up intermediate clutter (keeping the sorted BED for IGV)..."
rm -f "${demo_sam}" "${demo_bam}" "${demo_bam}.bai" "${spike_sam}" "${spike_bam}" "${merged_bam}"
rm -f "${input_bed_filename}" "${spike_bam_sorted%.bam}.bam" "${spike_bam_sorted}.bai"

echo "All done. Load:"
echo "  BAM: ${merged_bam_sorted}"
echo "  BED: ${output_bed_filename}"