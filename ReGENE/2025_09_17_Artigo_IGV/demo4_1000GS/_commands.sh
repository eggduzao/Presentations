#!/bin/bash
set -euo pipefail

# chr20:32,831,516-32,831,663 -> Insertion, Deletion and Prob. Het.

# FTP in MAC: sftp
lftp ftp.1000genomes.ebi.ac.uk

# Downloading from 1000genomes
wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/alignment/HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00100/alignment/HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam.bai

# Indexing and sorting bam
samtools index HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
samtools sort -o HG00096_chr20.bam HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
samtools index HG00096_chr20.bam

# Indexing and sorting bam
samtools index HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam
samtools sort -o HG00100_chr20.bam HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam
samtools index HG00100_chr20.bam

# Extracting and indexing vcf.gz - HG00096
bcftools view -s HG00096 ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -Oz -o HG00096.chr20.vcf.gz
tabix -p vcf HG00096.chr20.vcf.gz

# Extracting and indexing vcf.gz - HG00100
bcftools view -s HG00100 ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -Oz -o HG00100.chr20.vcf.gz
tabix -p vcf HG00100.chr20.vcf.gz

# Standardizing
mkdir -p _others/
mv HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam _others/
mv HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.bai _others/
mv HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam _others/
mv HG00100.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam.bai _others/

# Wrap all in bgzip
gzip -t  || echo "corrupt download"
gunzip -c ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | bgzip -c > ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.bgz

# Indexing the whole VCF
tabix -p vcf ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.bgz

# Standardizing
mkdir -p _allcvf/
mv ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz _allcvf/

# !! REMOVING EVERYTHING !!
rm -rf _allcvf/
rm -rf _others/
rm -rf HG00100_chr20.bam
rm -rf HG00100_chr20.bam.bai
rm -rf HG00100.chr20.vcf.gz
rm -rf ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.bgz
rm -rf ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.bgz.tbi
rm -rf HG00096_chr20.bam
rm -rf HG00096_chr20.bam.bai
rm -rf HG00096.chr20.vcf.gz
rm -rf HG00096.chr20.vcf.gz.tbi

