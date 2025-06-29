[![GitHub release](https://img.shields.io/github/v/release/Schahrjar/ExHap)](https://github.com/Schahrjar/ExHap/releases/latest)
[![last commit](https://img.shields.io/github/last-commit/Schahrjar/ExHap)](https://github.com/Schahrjar/ExHap/commits/main)

# ExHap
The tool retrievs homozygosity haplotypes from a multi-sample VCF file. Manily developed for exome variants of rare disorder patients to investigate if a common pathogenic mutation is a recurrent or founder variant. ExHap considers every signle variant in runs of homozygosity (ROH) of each sample to find any possible homozygous haplotype shared between two or more samples. Then reports the largest haplotype (and sub-haplotypes if any) in a BED file.

> [!NOTE]
> Feasibly find shared homozygous haplotypes in any cohort of genomic variants, visualise haplotypes, and calculate age of the most recent common acestor of any variant of interest.

## 🔧 Features
- Identifies all individual ROH per sample
- Detects homozygous haplotypes shared among samples
- Works on any genome assembly (e.g. hg38 or hg19)
- Requires Python

## 📦 Usage

## 🗂️ Inputs
### Mandatory
### Optional

## 📜 Citation

This tool is released prior to our manuscript submission to assist researchers. You may contact [Shahryar Alavi](https://schahrjar.github.io/) if you use ExHap in your publication.
