#!/bin/bash

# wget --content-disposition --trust-server-names \
#      -O K562_erythroleukemia_H3K4me3.bam \
#      "https://www.encodeproject.org/files/ENCFF197VMR/@@download/ENCFF197VMR.bam"

# Install samtools
# Windows: -
# Ubuntu: -
# MAC: brew install samtools

# Download Reads:
curl -L -o K562_erythroleukemia_H3K4me3.bam "https://www.encodeproject.org/files/ENCFF197VMR/@@download/ENCFF197VMR.bam"

# Download Signal over Control:
curl -L -o K562_erythroleukemia_H3K4me3_signal.bigWig "https://www.encodeproject.org/files/ENCFF814IYI/@@download/ENCFF814IYI.bigWig"

# Download Peaks Called:
curl -L -o K562_erythroleukemia_H3K4me3_peaks.bed.gz "https://www.encodeproject.org/files/ENCFF909PMV/@@download/ENCFF909PMV.bed.gz"
gunzip K562_erythroleukemia_H3K4me3_peaks.bed.gz

# Download Masked List:
curl -L -o hg38.chrom.sizes "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes"
curl -L -o hg38_Masked.bed.gz "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
gunzip hg38_Masked.bed.gz

# Sort bam file
samtools index K562_erythroleukemia_H3K4me3.bam
samtools sort -o K562_erythroleukemia_H3K4me3_sorted.bam K562_erythroleukemia_H3K4me3.bam

# Index with csi
samtools index -c K562_erythroleukemia_H3K4me3_sorted.bam

# Index with bai
samtools index K562_erythroleukemia_H3K4me3_sorted.bam

# !! DELETE EVERYTHING !!
rm -rf hg38_Masked_Sorted.bed
rm -rf K562_erythroleukemia_H3K4me3_sorted.bam
rm -rf K562_erythroleukemia_H3K4me3_sorted.bam.bai
rm -rf _others

