[![GitHub release](https://img.shields.io/github/v/release/Schahrjar/ExHap)](https://github.com/Schahrjar/ExHap/releases/latest)
[![last commit](https://img.shields.io/github/last-commit/Schahrjar/ExHap)](https://github.com/Schahrjar/ExHap/commits/main)

# ExHap
This tool identifies shared homozygosity haplotype blocks from a multi-sample VCF file. This is manily developed for investigating identity state of common pathogenic variants in rare disorder patients, if such variants are recurrent or founder. ExHap considers every signle variant in runs of homozygosity (ROH) of each sample to find any possible homozygous haplotypes shared between two or more samples. Then reports the largest haplotypes (and their sub-haplotypes if any) in a BED file.

ExHap has been tested on cohorts of exome data. It also visualises haplotypes and calculates age of the most recent common acestor of any variants of interest; if it is from a shared haplotype with the assumption that it is an identity by discent (IBD).

## 🔧 Features
- Detects homozygous haplotypes shared among samples
- Works on any genome assembly (e.g. hg38 or hg19)
- Requires Python

## 📦 Usage
Clone the repo

> [!WARNING]
> For WGS VCF data with hundreds of samples, it may need a large RAM.

## 🗂️ Inputs
### Mandatory
### Optional

## 📜 Citation

This tool is released prior to our manuscript submission to assist researchers. You may contact [Shahryar Alavi](https://schahrjar.github.io/) if you use ExHap in your publication.
