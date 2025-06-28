#!/usr/bin/env python3
import sys
import gzip
import argparse
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="ROH-DICE-lite: Group homozygous variants by sample pattern and proximity")
    parser.add_argument("-i", "--input", required=True, help="Input VCF (can be .vcf or .vcf.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    parser.add_argument("-d", "--maxdist", type=int, default=1_000_000, help="Max distance between variants in a group (default 1Mb)")
    return parser.parse_args()

def open_vcf(file):
    if file.endswith(".gz"):
        return gzip.open(file, "rt")
    else:
        return open(file, "r")

def extract_samples_and_genotypes(vcf_file):
    samples = []
    variants = []
    with open_vcf(vcf_file) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.strip().split('\t')
                samples = parts[9:]
                continue
            parts = line.strip().split('\t')
            chrom, pos, id_, ref, alt, qual, filt, info, fmt = parts[:9]
            genotypes = parts[9:]
            pos = int(pos)
            
            # Format fields are colon-separated, get GT field index:
            fmt_fields = fmt.split(":")
            try:
                gt_idx = fmt_fields.index("GT")
            except ValueError:
                # No GT field
                continue
            
            # Extract GT per sample
            gt_list = []
            for gt_field in genotypes:
                gt_data = gt_field.split(":")
                if len(gt_data) <= gt_idx:
                    gt_list.append("./.")  # missing
                else:
                    gt_list.append(gt_data[gt_idx])
            
            variants.append((chrom, pos, ref, alt, gt_list))
    return samples, variants

def is_homozygous(gt):
    # Accept only 0/0 or 1/1 (including phased and unphased)
    if gt in {"0/0", "0|0"}:
        return "HOM_REF"
    elif gt in {"1/1", "1|1"}:
        return "HOM_ALT"
    else:
        return None

def build_pattern(samples, gt_list):
    # For each sample, if homozygous, record sample:state, else ignore
    pattern_parts = []
    for s, gt in zip(samples, gt_list):
        state = is_homozygous(gt)
        if state is not None:
            pattern_parts.append(f"{s}:{state}")
    if pattern_parts:
        return "|".join(sorted(pattern_parts))
    else:
        return None

def group_variants(variants, maxdist):
    """
    Group variants by chrom and pattern, grouping adjacent variants
    of the same pattern within maxdist basepairs.
    """
    grouped = []
    # Data structure: chrom -> pattern -> list of (pos, ref, alt)
    chrom_pattern_pos = defaultdict(lambda: defaultdict(list))
    
    for chrom, pos, ref, alt, gt_list in variants:
        pattern = build_pattern(samples, gt_list)
        if pattern is None:
            continue
        chrom_pattern_pos[chrom][pattern].append((pos, ref, alt))
    
    for chrom in chrom_pattern_pos:
        for pattern in chrom_pattern_pos[chrom]:
            # Sort positions ascending
            pos_list = sorted(chrom_pattern_pos[chrom][pattern], key=lambda x: x[0])
            
            # Group nearby variants
            group_start = pos_list[0][0]
            group_end = pos_list[0][0]
            group_vars = [pos_list[0]]
            
            for i in range(1, len(pos_list)):
                curr_pos = pos_list[i][0]
                prev_pos = pos_list[i-1][0]
                if curr_pos - prev_pos <= maxdist:
                    group_end = curr_pos
                    group_vars.append(pos_list[i])
                else:
                    # save previous group
                    grouped.append((chrom, group_start, group_end, pattern, group_vars))
                    # start new group
                    group_start = curr_pos
                    group_end = curr_pos
                    group_vars = [pos_list[i]]
            # save last group
            grouped.append((chrom, group_start, group_end, pattern, group_vars))
    return grouped

def write_bed(groups, out_file):
    """
    BED fields:
    chrom (str)
    start (int, 0-based)
    end (int, 1-based)
    name (pattern string)
    score (number of variants in group)
    strand (.)
    details (semicolon separated variant positions 1-based)
    """
    with open(out_file, "w") as f:
        for chrom, start, end, pattern, vars_ in groups:
            # BED start is 0-based, so subtract 1 from start position
            bed_start = start - 1
            bed_end = end  # BED end is 1-based, non-inclusive but we keep 1-based here
            
            positions = [str(p[0]) for p in vars_]
            details = ";".join(positions)
            score = len(vars_)
            f.write(f"{chrom}\t{bed_start}\t{bed_end}\t{pattern}\t{score}\t.\t{details}\n")

if __name__ == "__main__":
    args = parse_args()
    samples, variants = extract_samples_and_genotypes(args.input)
    groups = group_variants(variants, args.maxdist)
    write_bed(groups, args.output)
    print(f"ROH-DICE-lite finished. Output written to {args.output}")
