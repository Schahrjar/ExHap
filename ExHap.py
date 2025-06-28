#!/usr/bin/env python3

import argparse
import gzip
from collections import defaultdict
from itertools import combinations
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Extract homozygous haplotypes from a VCF file and find shared haplotypes across samples.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF file (bgzipped and indexed).")
    parser.add_argument("-o", "--output", required=True, help="Output BED file for homozygous blocks.")
    parser.add_argument("-ho", "--haplotype_output", required=False, help="Output BED file for shared haplotypes.")
    parser.add_argument("-d", "--max_distance", type=int, default=10000, help="Max distance between variants to merge into a block.")
    parser.add_argument("-m", "--min_support", type=int, default=2, help="Minimum number of samples sharing a haplotype.")
    parser.add_argument("-l", "--min_loci", type=int, default=5, help="Minimum number of variant loci to define a shared haplotype.")
    return parser.parse_args()

def parse_vcf(input_file):
    genotypes = defaultdict(list)
    positions = []

    open_func = gzip.open if input_file.endswith(".gz") else open
    with open_func(input_file, 'rt') as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                samples = header[9:]
            elif line.startswith("#"):
                continue
            else:
                fields = line.strip().split("\t")
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                format_keys = fields[8].split(":")
                sample_fields = fields[9:]

                gt_index = format_keys.index("GT")
                block_genotypes = {}
                for sample, sample_data in zip(samples, sample_fields):
                    values = sample_data.split(":")
                    if len(values) <= gt_index:
                        continue
                    gt = values[gt_index].replace('|', '/')
                    if gt in {"0/0", "1/1"}:
                        block_genotypes[sample] = "HOM_REF" if gt == "0/0" else "HOM_ALT"
                genotypes[(chrom, pos)].append(block_genotypes)
                positions.append((chrom, pos))

    return samples, genotypes, positions

def cluster_blocks(positions, genotypes, max_distance):
    blocks = []
    block = []
    for i in range(len(positions)):
        chrom, pos = positions[i]
        if not block:
            block = [positions[i]]
        else:
            prev_chrom, prev_pos = block[-1]
            if chrom == prev_chrom and (pos - prev_pos <= max_distance):
                block.append(positions[i])
            else:
                blocks.append(block)
                block = [positions[i]]
    if block:
        blocks.append(block)
    return blocks

def write_bed(blocks, genotypes, samples, output_file):
    with open(output_file, 'w') as out:
        for block in blocks:
            chrom = block[0][0]
            start = block[0][1]
            end = block[-1][1]
            positions = [pos for _, pos in block]
            per_sample_gt = defaultdict(list)
            for pos in block:
                gt = genotypes[pos][0]
                for s in gt:
                    per_sample_gt[s].append(gt[s])
            genotype_summary = []
            for s in samples:
                if s in per_sample_gt:
                    if all(g == "HOM_REF" for g in per_sample_gt[s]):
                        genotype_summary.append(f"{s}:HOM_REF")
                    elif all(g == "HOM_ALT" for g in per_sample_gt[s]):
                        genotype_summary.append(f"{s}:HOM_ALT")
                    else:
                        genotype_summary.append(f"{s}:MIXED")
            allele_count = sum([1 for x in genotype_summary if "HOM_ALT" in x])
            out.write(f"{chrom}\t{start}\t{end}\t{'|'.join(genotype_summary)}\t{allele_count}\t.\t{';'.join(map(str, positions))}\n")

def extract_haplotypes(blocks, genotypes, samples, min_support, min_loci, haplo_out):
    haplo_dict = defaultdict(lambda: defaultdict(list))  # haplo_dict[sample][haplo_signature] = [(start, end, chrom, positions)]
    final_haplotypes = []

    for block in blocks:
        chrom = block[0][0]
        start = block[0][1]
        end = block[-1][1]
        positions = [pos for _, pos in block]

        for sample in samples:
            sig = []
            for pos in block:
                gt = genotypes[pos][0]
                sig.append(gt.get(sample, "NA"))
            if "NA" in sig:
                continue
            haplo_sig = ",".join(sig)
            haplo_dict[haplo_sig][tuple(sig)].append((start, end, chrom, positions, sample))

    for haplo_sig in haplo_dict:
        sample_map = defaultdict(list)
        for sig, entries in haplo_dict[haplo_sig].items():
            for start, end, chrom, positions, sample in entries:
                sample_map[(chrom, start, end, ",".join(map(str, positions)))].append(sample)

        for (chrom, start, end, pos_str), sample_list in sample_map.items():
            if len(sample_list) >= min_support and pos_str.count(',')+1 >= min_loci:
                final_haplotypes.append((chrom, start, end, sample_list, pos_str))

    with open(haplo_out, 'w') as out:
        for chrom, start, end, samples, pos_str in final_haplotypes:
            out.write(f"{chrom}\t{start}\t{end}\tHAPLOTYPE\t{';'.join(samples)}\t{len(samples)}\t{pos_str}\n")

def main():
    args = parse_args()
    samples, genotypes, positions = parse_vcf(args.input)
    blocks = cluster_blocks(positions, genotypes, args.max_distance)
    write_bed(blocks, genotypes, samples, args.output)

    if args.haplotype_output:
        extract_haplotypes(blocks, genotypes, samples, args.min_support, args.min_loci, args.haplotype_output)

if __name__ == "__main__":
    main()
