#!/bin/bash
set -euo pipefail

# Heart
# RNA-seq of human HEART tissue. The second technical replicate was performed in May 2015
# and was not part of the original Lin et al. 2014 publication (PMID:25413365).

# Isogenic replicate	1	1
# Technical replicate	1	2
# Summary	Homo sapiens female adult (30 years) heart tissue	Homo sapiens female adult (30 years) heart tissue
# Biosample	ENCBS227BAN	ENCBS227BAN
# Library	ENCLB797GKI	ENCLB171YLN
curl -L -o Heart_hg38_F30_TRep1.bam "https://www.encodeproject.org/files/ENCFF393JBM/@@download/ENCFF393JBM.bam" # Technical replicate 1
curl -L -o Heart_hg38_F30_TRep2.bam "https://www.encodeproject.org/files/ENCFF457BDV/@@download/ENCFF457BDV.bam" # Technical replicate 2

# Liver
# Homo sapiens liver tissue female child (4 years) and with
# nonobstructive coronary artery disease; liver tissue male adult (32 years)

# Anisogenic replicate	1	2
# Technical replicate	1	1
# Summary	Homo sapiens female child (4 years) liver tissue	Homo sapiens male adult (32 years) with nonobstructive coronary artery disease liver tissue
# Biosample	ENCBS401URL	ENCBS046RNA
# Library	ENCLB869UHR	ENCLB420PUQ
curl -L -o Liver_hg38_F04_TRep1.bam "https://www.encodeproject.org/files/ENCFF061GJD/@@download/ENCFF061GJD.bam" # Homo sapiens female child (4 years) liver tissue
curl -L -o Liver_hg38_M32_TRep1.bam "https://www.encodeproject.org/files/ENCFF699TZO/@@download/ENCFF699TZO.bam" # Homo sapiens male adult (32 years) with NCAD liver tissue

# Indexing all
samtools index Heart_hg38_F30_TRep1.bam
samtools index Heart_hg38_F30_TRep2.bam
samtools index Liver_hg38_F04_TRep1.bam
samtools index Heart_hg38_F30_TRep2.bam

# !! DELETING ALL !!
rm -rf Heart_hg38_F30_TRep1.bam
rm -rf Heart_hg38_F30_TRep1.bam.bai
rm -rf Heart_hg38_F30_TRep2.bam
rm -rf Heart_hg38_F30_TRep2.bam.bai
rm -rf Liver_hg38_F04_TRep1.bam
rm -rf Liver_hg38_F04_TRep1.bam.bai
rm -rf Liver_hg38_M32_TRep1.bam
rm -rf Liver_hg38_M32_TRep1.bam.bai

