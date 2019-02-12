#!/bin/bash

set -ex

REF=NC_016845.1
SAMPLE1=SRR3219271
SAMPLE2=SRR3219831

efetch -db nuccore -format fasta -id $REF > $REF.fasta

function run_snippy {
    local sample=$1
    fastq-dump --origfmt -I --split-files --gzip $sample
    trim_galore --paired --gzip "$sample"_*.fastq.gz
    snippy --outdir $sample --ref $REF.fasta --R1 "$sample"_1_val_1.fq.gz --R2 "$sample"_2_val_2.fq.gz
}

run_snippy "$SAMPLE1"
run_snippy "$SAMPLE2"

snippy-core --ref "$REF.fasta" "$SAMPLE1" "$SAMPLE2"

# get allele counts at sites where $SAMPLE1 != $SAMPLE2
bcftools view --min-ac 1:minor core.vcf -Oz -o diffs.vcf.gz
bcftools query -f '%CHROM\t%POS\n' diffs.vcf.gz > diffs.positions
bcftools mpileup --fasta-ref $REF.fasta -R diffs.positions $SAMPLE1/snps.bam $SAMPLE2/snps.bam -a FORMAT/AD -Oz -o repileup.vcf.gz

# merge snippy genotypes with mpileup allele depths
tabix diffs.vcf.gz
tabix repileup.vcf.gz
bcftools annotate -a repileup.vcf.gz -c FMT/AD diffs.vcf.gz > diffs_with_ad.vcf
