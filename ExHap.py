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

# Helper function to check if sub_list is a subsequence of main_list (order matters)
def is_subsequence(main_list, sub_list):
    if not sub_list:
        return True
    if not main_list:
        return False
    
    i = 0 # pointer for main_list
    j = 0 # pointer for sub_list
    while i < len(main_list) and j < len(sub_list):
        if main_list[i] == sub_list[j]:
            j += 1 
        i += 1 
    return j == len(sub_list)

def extract_haplotypes(blocks, genotypes, samples, min_support, min_loci, haplo_out):
    all_individual_homozygous_segments = []

    for block_of_positions in blocks: 
        current_chrom = block_of_positions[0][0] 

        for sample_name in samples: 
            current_segment_genotypes = []
            current_segment_variant_positions = []
            segment_start_genomic = -1 

            for pos_idx, (chrom_val, pos_genomic) in enumerate(block_of_positions):
                gt_at_pos = genotypes[(chrom_val, pos_genomic)][0].get(sample_name, None) 

                if gt_at_pos in {"HOM_REF", "HOM_ALT"}:
                    if not current_segment_genotypes: 
                        segment_start_genomic = pos_genomic
                    current_segment_genotypes.append(gt_at_pos)
                    current_segment_variant_positions.append(pos_genomic)
                else:
                    if current_segment_genotypes: 
                        if len(current_segment_genotypes) >= min_loci:
                            segment_end_genomic = current_segment_variant_positions[-1]
                            all_individual_homozygous_segments.append(
                                (current_chrom, segment_start_genomic, segment_end_genomic,
                                 tuple(current_segment_genotypes), tuple(current_segment_variant_positions), sample_name)
                            )
                        current_segment_genotypes = []
                        current_segment_variant_positions = []
                        segment_start_genomic = -1 
            
            if current_segment_genotypes:
                if len(current_segment_genotypes) >= min_loci:
                    segment_end_genomic = current_segment_variant_positions[-1]
                    all_individual_homozygous_segments.append(
                        (current_chrom, segment_start_genomic, segment_end_genomic,
                         tuple(current_segment_genotypes), tuple(current_segment_variant_positions), sample_name)
                    )

    all_valid_homozygous_subsegments = []

    for seg_chrom, _, _, seg_genotypes_orig, seg_variants_orig, sample_name in all_individual_homozygous_segments:
        num_loci_in_orig_segment = len(seg_genotypes_orig)
        
        for i in range(num_loci_in_orig_segment):
            for j in range(i, num_loci_in_orig_segment):
                sub_genotypes = seg_genotypes_orig[i : j+1]
                sub_variants = seg_variants_orig[i : j+1]

                if len(sub_genotypes) >= min_loci:
                    sub_start_genomic = sub_variants[0]
                    sub_end_genomic = sub_variants[-1]
                    
                    all_valid_homozygous_subsegments.append(
                        (seg_chrom, sub_start_genomic, sub_end_genomic,
                         sub_genotypes, sub_variants, sample_name)
                    )

    shared_haplotype_candidate = defaultdict(list)
    for sub_chrom, sub_start, sub_end, sub_genotypes, sub_variants, sample_name in all_valid_homozygous_subsegments:
        key = (sub_chrom, sub_start, sub_end, sub_genotypes, sub_variants)
        shared_haplotype_candidate[key].append(sample_name)

    # Convert to a list of dicts for easier sorting and processing
    temp_final_haplotypes = []
    for (chrom, start, end, genotypes_tuple, variants_tuple), sample_list in shared_haplotype_candidate.items():
        if len(sample_list) >= min_support:
            temp_final_haplotypes.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'genotypes': genotypes_tuple, # Keep genotype pattern for comparison if needed in future
                'positions': variants_tuple,
                'samples': frozenset(sample_list), # Use frozenset for efficient comparison
                'original_sample_list': sample_list, # Keep original list for output
                'original_pos_str': ",".join(map(str, variants_tuple)) # Keep original string for output
            })
    
    # MODIFICATION START: Filter redundant sub-segments (Point 3)
    # Sort haplotypes: longest first, then by chromosome, then by start position
    temp_final_haplotypes.sort(key=lambda x: (x['chrom'], -len(x['positions']), x['start']))

    filtered_haplotypes = []
    redundant_indices = set() # Store indices of haplotypes that are found to be redundant

    for i in range(len(temp_final_haplotypes)):
        if i in redundant_indices:
            continue # This haplotype has already been marked as redundant

        h1 = temp_final_haplotypes[i] # This is a candidate for a non-redundant (often longer) haplotype
        
        # Add h1 to the filtered list
        filtered_haplotypes.append(h1)

        # Compare h1 with all subsequent haplotypes (h2)
        for j in range(i + 1, len(temp_final_haplotypes)):
            if j in redundant_indices:
                continue # This haplotype has already been marked as redundant

            h2 = temp_final_haplotypes[j] # This is a candidate for being a redundant (shorter) haplotype

            # Check conditions for redundancy:
            # 1. Same chromosome
            # 2. Exact same set of samples
            # 3. h2's positions are a subsequence of h1's positions (meaning h2 is a sub-segment of h1)
            if (h1['chrom'] == h2['chrom'] and
                h1['samples'] == h2['samples'] and
                is_subsequence(h1['positions'], h2['positions'])):
                
                redundant_indices.add(j) # Mark h2 as redundant
    # MODIFICATION END

    haplotype_number = 0 
    with open(haplo_out, 'w') as out:
        for hap_data in filtered_haplotypes:
            haplotype_number += 1 
            haplotype_name = f"HAPLOTYPE{haplotype_number}" 
            out.write(f"{hap_data['chrom']}\t{hap_data['start']}\t{hap_data['end']}\t"
                      f"{haplotype_name}\t{len(hap_data['original_sample_list'])}\t"
                      f"{';'.join(hap_data['original_sample_list'])}\t{hap_data['original_pos_str']}\n")

def main():
    args = parse_args()
    samples, genotypes, positions = parse_vcf(args.input)
    blocks = cluster_blocks(positions, genotypes, args.max_distance)
    write_bed(blocks, genotypes, samples, args.output)

    if args.haplotype_output:
        extract_haplotypes(blocks, genotypes, samples, args.min_support, args.min_loci, args.haplotype_output)

if __name__ == "__main__":
    main()
