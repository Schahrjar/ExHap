#!/usr/bin/env python3
import argparse
import gzip
import sys
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="ExHap: detect shared homozygous haplotypes.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF (.vcf or .vcf.gz)")
    parser.add_argument("-o", "--output", required=True, help="BED file of homozygous blocks")
    parser.add_argument("-d", "--maxdist", type=int, default=10000, help="Max distance between variants in a block")
    parser.add_argument("-ho", "--haplotype_output", help="BED file for shared haplotypes")
    parser.add_argument("--min_hap_length", type=int, default=0, help="Minimum haplotype length (bp) to output")
    parser.add_argument("--min_support", type=int, default=2, help="Minimum samples sharing a haplotype")
    parser.add_argument("--min_loci", type=int, default=2, help="Minimum variant count defining a haplotype")
    return parser.parse_args()

def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def extract_samples_and_genotypes(vcf_file):
    samples = []
    variants = []
    with open_vcf(vcf_file) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                continue
            fields = line.strip().split("\t")
            chrom, pos, _, ref, alt, *_ , fmt = fields[:9]
            pos = int(pos)
            gt_index = fmt.split(":").index("GT")
            gt_list = [s.split(":")[gt_index] if len(s.split(":")) > gt_index else "./." for s in fields[9:]]
            variants.append((chrom, pos, gt_list))
    return samples, variants

def is_homozygous(gt):
    if gt in ("0/0", "0|0"):
        return "HOM_REF"
    if gt in ("1/1", "1|1"):
        return "HOM_ALT"
    return None

def build_pattern(samples, gt_list):
    parts = []
    for s, gt in zip(samples, gt_list):
        h = is_homozygous(gt)
        if h:
            parts.append(f"{s}:{h}")
    return "|".join(sorted(parts)) if parts else None

def group_variants(variants, maxdist, samples):
    chrom_pattern = defaultdict(lambda: defaultdict(list))
    for chrom, pos, gt_list in variants:
        pat = build_pattern(samples, gt_list)
        if not pat:
            continue
        chrom_pattern[chrom][pat].append(pos)

    groups = []
    for chrom in chrom_pattern:
        for pat in chrom_pattern[chrom]:
            pos_list = sorted(chrom_pattern[chrom][pat])
            start = pos_list[0]
            end = start
            temp = [start]
            for p in pos_list[1:]:
                if p - end <= maxdist:
                    end = p
                    temp.append(p)
                else:
                    groups.append((chrom, start, end, pat, temp))
                    start = p
                    end = p
                    temp = [p]
            groups.append((chrom, start, end, pat, temp))
    return groups

def write_bed(groups, out):
    for chrom, start, end, pat, pos_list in groups:
        score = len(pos_list)
        out.write(f"{chrom}\t{start-1}\t{end}\t{pat}\t{score}\t.\t{';'.join(map(str,pos_list))}\n")

def extract_shared_haplotypes(groups, samples, args):
    haps = []
    hid = 1
    for i, (chrom, s0, e0, pat0, pos0) in enumerate(groups):
        samples0 = set(x.split(":")[0] for x in pat0.split("|"))
        gt0 = pat0.split("|")[0].split(":")[1]  # homozygous type
        if len(samples0) < args.min_support or len(pos0) < args.min_loci: 
            continue

        start, end, positions = s0, e0, pos0.copy()
        shared = samples0.copy()
        for _, s1, e1, pat1, pos1 in groups[i+1:]:
            if _ != chrom or s1 - end > args.maxdist:
                break
            samples1 = set(x.split(":")[0] for x in pat1.split("|"))
            merged = shared & samples1
            if len(merged) < args.min_support:
                break
            shared = merged
            end = e1
            positions.extend(pos1)

        length = end - start + 1
        if shared and len(shared) >= args.min_support and len(positions) >= args.min_loci and length >= args.min_hap_length:
            haps.append((chrom, start, end, f"HAPLOTYPE{hid}", len(shared), length, ",".join(sorted(shared))))
            hid += 1
    return haps

def main():
    args = parse_args()
    samples, variants = extract_samples_and_genotypes(args.input)
    groups = group_variants(variants, args.maxdist, samples)

    with open(args.output, "w") as bed_out:
        write_bed(groups, bed_out)
    print(f"Homozygous blocks written to {args.output}")

    if args.haplotype_output:
        haps = extract_shared_haplotypes(groups, samples, args)
        with open(args.haplotype_output, "w") as hap_out:
            for rec in haps:
                hap_out.write("\t".join(map(str, rec))+"\n")
        print(f"Shared haplotypes written to {args.haplotype_output} ({len(haps)} entries)")

if __name__ == "__main__":
    main()
