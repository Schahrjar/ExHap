import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys
import argparse
from scipy.stats import gamma # For MRCA calculation
import matplotlib.ticker as mticker # For integer colorbar ticks
import matplotlib.patches as mpatches # For legend handles

# --- MRCA Calculation Function ---
def calculate_mrca_age(l_lengths_mbp, r_lengths_mbp, confidence_coefficient=0.95, generation_time_years=20):
    """
    Calculates MRCA age based on shared homozygous haplotype lengths.
    Translates logic from the provided R script.
    Lengths are expected in Mbp, which are treated as cM (centiMorgans).
    Internally, these are converted to Morgans as per the R script.
    
    Args:
        l_lengths_mbp (list): Genetic distances (in Mbp) from variant to left end of haplotype.
        r_lengths_mbp (list): Genetic distances (in Mbp) from right end of haplotype.
        confidence_coefficient (float): Desired confidence level (e.g., 0.95 for 95% CI).
        generation_time_years (int/float): Number of years per generation.

    Returns:
        tuple: (independent_age_generations, (independent_l_ci_gen, independent_u_ci_gen),
                independent_age_years, (independent_l_ci_years, independent_u_ci_years),
                correlated_age_generations, (correlated_l_ci_gen, correlated_u_ci_gen),
                correlated_age_years, (correlated_l_ci_years, correlated_u_ci_years))
               Returns (None, None, None, None, None, None, None, None) if not enough data or calculation fails.
    """
    l_lengths_cm = np.array(l_lengths_mbp) # Treat Mbp as cM
    r_lengths_cm = np.array(r_lengths_mbp)

    if len(l_lengths_cm) != len(r_lengths_cm) or len(l_lengths_cm) == 0:
        print("MRCA Warning: l_lengths and r_lengths must have the same non-zero number of elements for MRCA calculation.")
        return None, None, None, None, None, None, None, None

    # Convert cM to Morgans as per R script: l.lengths = (1 / 100) * l.lengths
    l_lengths = l_lengths_cm / 100 
    r_lengths = r_lengths_cm / 100

    n = len(l_lengths)
    cc = confidence_coefficient

    # As in the provided R script example, chance.sharing.correction is FALSE
    cs_correction = 0 
    # In the R script, i.cs.correction is 0 if n < 10, otherwise cs.correction.
    # Given cs_correction = 0 (from example), i.cs.correction will always be 0.
    i_cs_correction = 0 
    
    # Initialize all return values with None or "N/A"
    independent_tau_hat, independent_l, independent_u = "N/A", "N/A", "N/A"
    independent_tau_hat_years, independent_l_years, independent_u_years = "N/A", "N/A", "N/A"
    correlated_tau_hat, correlated_l, correlated_u = "N/A", "N/A", "N/A"
    correlated_tau_hat_years, correlated_l_years, correlated_u_years = "N/A", "N/A", "N/A"

    try:
        # ===== Age estimation - Assuming an 'independent' genealogy =====
        length_correction_independent = (np.sum(l_lengths) + np.sum(r_lengths) - 2 * (n - 1) * i_cs_correction) / (2 * n)
        sum_lengths_independent = np.sum(l_lengths) + np.sum(r_lengths) + 2 * length_correction_independent - 2 * (n - 1) * i_cs_correction
        
        if sum_lengths_independent > 0:
            b_c_independent = (2 * n - 1) / (2 * n)
            independent_tau_hat = (b_c_independent * 2 * n) / sum_lengths_independent
            
            g_l_independent = gamma.ppf(((1 - cc) / 2), a=2 * n, scale=1 / (2 * n * b_c_independent))
            g_u_independent = gamma.ppf((cc + (1 - cc) / 2), a=2 * n, scale=1 / (2 * n * b_c_independent))
            independent_l = g_l_independent * independent_tau_hat
            independent_u = g_u_independent * independent_tau_hat
    except Exception as e:
        print(f"MRCA Warning: Error in independent genealogy calculation: {e}")

    try:
        # ===== Age estimation - Assuming a 'correlated' genealogy =====
        length_correction_correlated = (np.sum(l_lengths) + np.sum(r_lengths) - 2 * (n - 1) * cs_correction) / (2 * n)

        l_lengths_corr = np.copy(l_lengths)
        r_lengths_corr = np.copy(r_lengths)

        if n > 0: # Ensure there are elements to find max
            longest_l_idx = np.argmax(l_lengths_corr)
            l_lengths_corr[longest_l_idx] += length_correction_correlated + cs_correction
            longest_r_idx = np.argmax(r_lengths_corr)
            r_lengths_corr[longest_r_idx] += length_correction_correlated + cs_correction

        lengths = l_lengths_corr + r_lengths_corr
        lengths = lengths - 2 * cs_correction 
        
        rho_hat = 0
        n_star = n # Default if n=1 or issues
        if n > 1: # Only calculate rho_hat for n > 1
            mean_lengths = np.mean(lengths)
            var_lengths = np.var(lengths)

            if mean_lengths == 0 and var_lengths == 0: # All lengths are zero
                 rho_hat = 1.0 # Perfect correlation
            elif var_lengths == 0: # No variance in lengths, but mean is non-zero
                rho_hat = 1.0 # Perfect correlation
            elif mean_lengths == 0: # Mean is zero, but variance is not (shouldn't happen with positive lengths)
                 rho_hat = 0.0 # Cannot compute, assume no correlation
            else:
                numerator = (n * (mean_lengths) ** 2 - var_lengths * (1 + 2 * n))
                denominator = (n * (mean_lengths) ** 2 + var_lengths * (n - 1))
                if denominator == 0:
                    rho_hat = 0.0 # Cannot compute, assume no correlation
                else:
                    rho_hat = numerator / denominator

            n_star = n / (1 + (n - 1) * rho_hat)
            # Apply N* capping as in R script, and ensure positive
            if n_star > n:
                n_star = n
            elif n_star <= 0: # Must be positive for gamma distribution
                n_star = 1 # Minimum effective sample size
            
            # Additional correction for n_star based on rho_hat in R script
            if rho_hat < -2 / (n - 1):
                n_star = n / (1 + (n - 1) * abs(rho_hat))
            elif -2/(n-1) <= rho_hat and rho_hat < -1/(n-1):
                n_star = n
            
            if n_star <= 0: # Double check, should be caught above but for safety
                n_star = 1

        if np.sum(lengths) > 0 and n_star > 0:
            b_c_correlated = (2 * n_star - 1) / (2 * n_star)
            correlated_tau_hat = (b_c_correlated * 2 * n) / np.sum(lengths) # R script uses 2*n here for num, not 2*n_star
            
            g_l_correlated = gamma.ppf(((1 - cc) / 2), a=2 * n_star, scale=1 / (2 * n_star * b_c_correlated))
            g_u_correlated = gamma.ppf((cc + (1 - cc) / 2), a=2 * n_star, scale=1 / (2 * n_star * b_c_correlated))
            correlated_l = g_l_correlated * correlated_tau_hat
            correlated_u = g_u_correlated * correlated_tau_hat
    except Exception as e:
        print(f"MRCA Warning: Error in correlated genealogy calculation: {e}")

    # Ensure finite numbers for output and handle inf/nan gracefully for rounding
    independent_tau_hat = round(independent_tau_hat, 1) if np.isfinite(independent_tau_hat) else "N/A"
    independent_l = round(independent_l, 1) if np.isfinite(independent_l) else "N/A"
    independent_u = round(independent_u, 1) if np.isfinite(independent_u) else "N/A"

    correlated_tau_hat = round(correlated_tau_hat, 1) if np.isfinite(correlated_tau_hat) else "N/A"
    correlated_l = round(correlated_l, 1) if np.isfinite(correlated_l) else "N/A"
    correlated_u = round(correlated_u, 1) if np.isfinite(correlated_u) else "N/A"

    # Calculate years
    if isinstance(independent_tau_hat, (int, float)):
        independent_tau_hat_years = round(independent_tau_hat * generation_time_years, 0)
        independent_l_years = round(independent_l * generation_time_years, 0) if isinstance(independent_l, (int, float)) else "N/A"
        independent_u_years = round(independent_u * generation_time_years, 0) if isinstance(independent_u, (int, float)) else "N/A"
    
    if isinstance(correlated_tau_hat, (int, float)):
        correlated_tau_hat_years = round(correlated_tau_hat * generation_time_years, 0)
        correlated_l_years = round(correlated_l * generation_time_years, 0) if isinstance(correlated_l, (int, float)) else "N/A"
        correlated_u_years = round(correlated_u * generation_time_years, 0) if isinstance(correlated_u, (int, float)) else "N/A"

    return (independent_tau_hat, (independent_l, independent_u), independent_tau_hat_years, (independent_l_years, independent_u_years),
            correlated_tau_hat, (correlated_l, correlated_u), correlated_tau_hat_years, (correlated_l_years, correlated_u_years))

def plot_shared_haplotypes(haplotype_output_file, target_chrom=None, genomic_range=None, mrca_variant_str=None, gene_name=None, color_by='samples'):
    """
    Plots shared homozygous haplotypes on a genomic scale and saves to file.
    Haplotypes on the same chromosome will be stacked vertically.
    Color indicates the number of samples sharing the haplotype.
    Includes optional range filtering and MRCA age calculation/display.
    Sub-segments are visually separated by sorting.
    """
    try:
        df = pd.read_csv(haplotype_output_file, sep='\t', header=None, 
                         names=['chrom', 'start', 'end', 'Haplotype_ID', 'Num_Samples', 'Sample_Names', 'Variant_Loci_Positions'])
    except FileNotFoundError:
        print(f"Error: Input file '{haplotype_output_file}' not found.")
        return
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    # --- Determine the plotting region based on arguments ---
    plot_chrom = target_chrom
    plot_start = None
    plot_end = None

    if genomic_range:
        try:
            parts = genomic_range.split(':')
            plot_chrom = parts[0]
            start_end_parts = parts[1].split('-')
            plot_start = int(start_end_parts[0])
            plot_end = int(start_end_parts[1])
            if plot_start >= plot_end:
                raise ValueError("Start position must be less than end position.")
            print(f"DEBUG: Plotting within specified range: {plot_chrom}:{plot_start:,}-{plot_end:,}")
        except (IndexError, ValueError):
            print(f"Error: Invalid genomic range format '{genomic_range}'. Expected 'chr:start-end'. Plotting will proceed without range filtering.")
            plot_start, plot_end = None, None # Disable range filtering

    # Filter DataFrame by chromosome
    if plot_chrom:
        df = df[df['chrom'] == plot_chrom]
        if df.empty:
            print(f"No haplotypes found for chromosome {plot_chrom} in '{haplotype_output_file}'. No plot will be generated.")
            return

    # Further filter by genomic start/end if a range was provided
    if plot_start is not None and plot_end is not None:
        # Haplotypes must overlap the range
        df = df[(df['start'] < plot_end) & (df['end'] > plot_start)] 
        if df.empty:
            print(f"No haplotypes found within range {plot_chrom}:{plot_start:,}-{plot_end:,} after filtering. No plot will be generated.")
            return

    # Sort for better visualization (by chromosome, then number of samples desc, then start position, then end)
    # This sorting order helps to visualize sub-segments: larger groups (more samples) and longer segments
    # will generally appear at lower y-levels, while smaller, nested segments will appear higher.
    df = df.sort_values(by=['chrom', 'Num_Samples', 'start', 'end'], ascending=[True, False, True, True])

    # Determine unique chromosomes for plotting (should be just one if chrom filtered)
    unique_chroms = df['chrom'].unique()
    
    if df.empty:
        print("No data to plot after all filtering. No plot will be generated.")
        return

    # --- MRCA Calculation Setup ---
    mrca_variant_chrom, mrca_variant_pos = None, None
    mrca_result_str = ""
    GEN_TIME_YEARS = 20 # Constant for generation time

    if mrca_variant_str:
        try:
            mrca_parts = mrca_variant_str.split(':')
            mrca_variant_chrom = mrca_parts[0]
            mrca_variant_pos = int(mrca_parts[1])
            print(f"DEBUG: MRCA variant specified for calculation: {mrca_variant_chrom}:{mrca_variant_pos:,}")

            if plot_chrom and mrca_variant_chrom != plot_chrom:
                print(f"MRCA Warning: MRCA variant chromosome '{mrca_variant_chrom}' does not match plotting chromosome '{plot_chrom}'. MRCA will not be calculated for this plot.")
                mrca_variant_str = None # Disable MRCA if chromosomes don't match
            elif plot_start is not None and plot_end is not None and \
                 not (plot_start <= mrca_variant_pos <= plot_end):
                print(f"MRCA Warning: MRCA variant {mrca_variant_chrom}:{mrca_variant_pos:,} is outside the plotting range {plot_start:,}-{plot_end:,}. MRCA will not be calculated for this plot.")
                mrca_variant_str = None # Disable MRCA if outside plot range
            else:
                mrca_l_lengths = []
                mrca_r_lengths = []
                
                # Iterate over the *already filtered* DataFrame 'df'
                for _, row in df.iterrows(): 
                    # Ensure haplotype segment contains the MRCA variant
                    if row['start'] <= mrca_variant_pos <= row['end']:
                        mrca_l_lengths.append((mrca_variant_pos - row['start']) / 1_000_000) # Convert bp to Mbp
                        mrca_r_lengths.append((row['end'] - mrca_variant_pos) / 1_000_000) # Convert bp to Mbp

                if mrca_l_lengths:
                    print(f"DEBUG: Found {len(mrca_l_lengths)} shared homozygous haplotype segments spanning {mrca_variant_str} for MRCA calculation.")
                    
                    indep_gen, indep_gen_ci, indep_year, indep_year_ci, \
                    corr_gen, corr_gen_ci, corr_year, corr_year_ci = \
                        calculate_mrca_age(mrca_l_lengths, mrca_r_lengths, generation_time_years=GEN_TIME_YEARS)

                    if indep_gen is not None:
                        mrca_result_str = (f"MRCA (Independent):\n"
                                           f"  {indep_gen} gens ({indep_gen_ci[0]}-{indep_gen_ci[1]} CI)\n"
                                           f"  {indep_year} years ({indep_year_ci[0]}-{indep_year_ci[1]} CI)\n"
                                           f"MRCA (Correlated):\n"
                                           f"  {corr_gen} gens ({corr_gen_ci[0]}-{corr_gen_ci[1]} CI)\n"
                                           f"  {corr_year} years ({corr_year_ci[0]}-{corr_year_ci[1]} CI)")
                        print(f"MRCA Calculation Result:\n{mrca_result_str}")
                    else:
                        print("MRCA Warning: MRCA calculation returned no valid results. Check debug messages above.")
                else:
                    print(f"MRCA Warning: No shared homozygous haplotype segments found spanning {mrca_variant_str} in the filtered data. Cannot calculate MRCA.")

        except ValueError:
            print(f"MRCA Error: Invalid MRCA variant format '{mrca_variant_str}'. Expected 'chrom:pos'. Skipping MRCA calculation.")
            mrca_variant_str = None # Disable MRCA calculation

    # --- Coloring Setup ---
    cmap = None
    norm = None
    cbar_label = ''
    show_colorbar = False
    show_legend = False

    if color_by == 'samples':
        min_samples = df['Num_Samples'].min()
        max_samples = df['Num_Samples'].max()
        if min_samples == max_samples:
            norm = plt.Normalize(vmin=min_samples, vmax=min_samples + 1)
        else:
            norm = plt.Normalize(vmin=min_samples, vmax=max_samples)
        cmap = cm.viridis # For continuous scale
        cbar_label = 'Number of Samples'
        show_colorbar = True
        show_legend = False

    elif color_by == 'haplotype_id':
        unique_haplotype_ids = sorted(df['Haplotype_ID'].unique())
        num_unique_haplotypes = len(unique_haplotype_ids)
        
        # Create a mapping from Haplotype_ID to an integer index
        hap_id_to_idx = {hap_id: i for i, hap_id in enumerate(unique_haplotype_ids)}
        
        # Choose a discrete colormap
        if num_unique_haplotypes <= 20:
            cmap = plt.get_cmap('tab20')
        else: # For >20 unique IDs, use a colormap that cycles colors
            cmap = plt.get_cmap('hsv') # Hsv is good for distinct categorical colors
        
        # Normalize to the range [0, num_unique_haplotypes - 1] for the cmap
        norm = plt.Normalize(vmin=0, vmax=num_unique_haplotypes - 1) 

        cbar_label = 'Haplotype ID' # Will be used if a colorbar is shown (but here we'll use a legend)
        show_colorbar = False 
        show_legend = True

    # --- Plotting Setup ---
    # Calculate initial y_level map to determine max y_level for plot height
    y_level_map_temp = {}
    current_y_level_temp = 0
    for idx, row in df.iterrows():
        if row['Haplotype_ID'] not in y_level_map_temp:
            y_level_map_temp[row['Haplotype_ID']] = current_y_level_temp
            current_y_level_temp += 1

    fig_height = max(5, current_y_level_temp * 0.4 + 2) # Adjusted dynamic height
    fig, axes = plt.subplots(len(unique_chroms), 1, figsize=(12, fig_height), 
                             sharex=False, squeeze=False) 
    axes = axes.flatten() 

    for i, chrom in enumerate(unique_chroms): 
        ax = axes[i]
        chrom_df = df[df['chrom'] == chrom]

        y_level_map = {} # Reset for each chromosome
        current_y_level = 0

        # Plot haplotypes
        for idx, row in chrom_df.iterrows():
            if row['Haplotype_ID'] not in y_level_map:
                y_level_map[row['Haplotype_ID']] = current_y_level
                current_y_level += 1
            y_level = y_level_map[row['Haplotype_ID']]
            
            # Determine color based on chosen method
            if color_by == 'samples':
                color = cmap(norm(row['Num_Samples']))
            elif color_by == 'haplotype_id':
                hap_idx = hap_id_to_idx[row['Haplotype_ID']]
                color = cmap(norm(hap_idx))
            
            ax.barh(y_level, width=row['end'] - row['start'], left=row['start'], 
                    height=0.8, color=color, edgecolor='black', linewidth=0.5)
            
        ax.set_title(f"Chromosome {chrom}")
        ax.set_ylabel("Haplotype Band")
        ax.set_xlabel("Genomic Position (bp)")
        ax.set_yticks([]) 
        if y_level_map:
            ax.set_ylim(-1, max(y_level_map.values()) + 1)
        else:
            ax.set_ylim(-1, 1)

        ax.ticklabel_format(style='plain', axis='x', useOffset=False) 
        ax.get_xaxis().get_major_formatter().set_scientific(False)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

        # Set X-limits based on the range if provided, otherwise auto
        if plot_start is not None and plot_end is not None and chrom == plot_chrom:
            ax.set_xlim(plot_start, plot_end)
        else:
            min_x = chrom_df['start'].min()
            max_x = chrom_df['end'].max()
            ax.set_xlim(min_x, max_x)

        # Plot vertical line at MRCA variant position if applicable to this subplot
        if mrca_variant_str and chrom == mrca_variant_chrom and \
           ( (plot_start is None and plot_end is None) or \
             (plot_start <= mrca_variant_pos <= plot_end) ):
            
            ax.axvline(x=mrca_variant_pos, color='red', linestyle='--', linewidth=1.5, zorder=3)
            # Add the actual variant string (e.g., 'chr4:185,415,857') as the label
            ax.text(mrca_variant_pos, ax.get_ylim()[1] * 0.9, mrca_variant_str, 
                    rotation=90, va='top', ha='right', color='red', fontsize=8, zorder=4)


    # --- Color Bar / Legend Section ---
    rect_val_right_adj = 0.82 # Default for colorbar
    if show_colorbar:
        fig.subplots_adjust(right=rect_val_right_adj)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7]) 
        cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax)
        cbar.set_label(cbar_label) 

        # Set informative integer ticks for the color bar
        if min_samples <= max_samples:
            unique_sample_counts = sorted(df['Num_Samples'].unique())
            if len(unique_sample_counts) <= 15: # If few unique values, show all
                cbar.set_ticks(unique_sample_counts)
                cbar.set_ticklabels([str(int(t)) for t in unique_sample_counts])
            else: # Otherwise, use MaxNLocator to select a reasonable number
                locator = mticker.MaxNLocator(nbins=10, integer=True, prune='both')
                cbar.locator = locator
                cbar.update_ticks()
                chosen_ticks = cbar.get_ticks()
                # Ensure labels are integers for chosen ticks
                if len(chosen_ticks) > 0:
                    cbar.set_ticklabels([str(int(t)) for t in chosen_ticks])
                elif min_samples == max_samples: # Edge case: only one value
                    cbar.set_ticks([min_samples])
                    cbar.set_ticklabels([str(int(min_samples))])
    
    if show_legend:
        # Create a proxy artist for each unique Haplotype_ID for the legend
        legend_elements = []
        for hap_id in unique_haplotype_ids:
            hap_idx = hap_id_to_idx[hap_id]
            color = cmap(norm(hap_idx))
            legend_elements.append(mpatches.Rectangle((0, 0), 1, 1, fc=color, ec="black", 
                                                      label=str(hap_id))) # Use Haplotype ID string directly

        # Adjust right margin to make space for the legend
        rect_val_right_adj = 0.75 # More space for legend than for colorbar
        fig.subplots_adjust(right=rect_val_right_adj)
        
        # Position the legend outside the main plot area to the right
        fig.legend(handles=legend_elements, title="Haplotype IDs", 
                   loc='center right', bbox_to_anchor=(0.99, 0.5), 
                   fontsize='small', bbox_transform=fig.transFigure)
        
        if num_unique_haplotypes > 20: # Warn if many IDs, as colors will repeat for tab20
            print("Warning: More than 20 unique Haplotype IDs. Colors in legend may repeat if 'tab20' is used, or vary continuously if 'hsv' is used.")


    # Adjust tight_layout
    # Ensure left margin respects potential MRCA text
    rect_left_adj = 0.05 if mrca_result_str else 0.0
    
    plt.tight_layout(rect=[rect_left_adj, 0, rect_val_right_adj, 1]) 
    
    # Add MRCA results as a text box on the figure
    if mrca_result_str:
        # Position relative to figure (0 to 1)
        fig.text(0.01, 0.99, "MRCA Age Estimation:\n" + mrca_result_str, 
                 verticalalignment='top', horizontalalignment='left', 
                 transform=fig.transFigure,
                 fontsize=9, bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="gray", lw=1, alpha=0.8))

    # --- Set plot title based on gene_name ---
    plot_title = "Shared Homozygous Haplotypes"
    if gene_name:
        plot_title = f"Shared Homozygous Haplotypes around {gene_name}"

    plt.suptitle(plot_title, y=1.02, fontsize=16) 

    # --- Determine output filename based on arguments ---
    filename_parts = []
    if gene_name:
        filename_parts.append(gene_name.replace(" ", "_")) # Use gene name at the start
    else:
        filename_parts.append("shared_haplotypes_plot")

    if genomic_range:
        filename_parts.append(genomic_range.replace(':', '_').replace('-', '_'))
    elif plot_chrom: 
        filename_parts.append(plot_chrom)
    if mrca_variant_str:
        filename_parts.append(f"MRCA_{mrca_variant_str.replace(':', '_')}")
    
    # Add coloring method to filename for distinction
    filename_parts.append(f"colorby_{color_by}")

    output_filename = "_".join(filename_parts) + ".png"
        
    plt.savefig(output_filename, dpi=300, bbox_inches='tight') 
    print(f"Plot saved to {output_filename}")


# --- Command-line argument parsing and function call ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot shared homozygous haplotypes from a BED-like file, with optional genomic range filtering and MRCA age calculation.",
        formatter_class=argparse.RawTextHelpFormatter 
    )

    parser.add_argument("haplotype_output_file", help="Path to the shared haplotype output file (from -ho argument).")
    
    parser.add_argument("--range", 
                        help="[Optional] Genomic range to plot, in 'chrom:start-end' format (e.g., 'chr4:100000-200000'). "
                             "If provided, filters the plot to this specific region.")
    
    parser.add_argument("--mrca_variant", 
                        help="[Optional] A specific variant locus (e.g., 'chr4:1234567') to calculate and display MRCA age for. "
                             "MRCA is calculated based on shared haplotypes spanning this variant. "
                             "Requires 'scipy' library for calculation.")

    parser.add_argument("--gene_name", 
                        help="[Optional] Name of the gene of interest to include in the plot title and output filename. "
                             "E.g., 'BRCA1'.")

    parser.add_argument("--color_by", 
                        choices=['samples', 'haplotype_id'],
                        default='samples',
                        help="[Optional] Method to color the haplotypes:\n"
                             "- 'samples': Color by the number of samples sharing the haplotype (default).\n"
                             "- 'haplotype_id': Color by the unique Haplotype_ID, representing its variant pattern.")

    # Positional target_chromosome for backward compatibility
    parser.add_argument("target_chromosome", nargs='?', 
                        help="[Optional] A single chromosome (e.g., 'chr4') to plot. "
                             "If --range is provided, this argument is ignored for filtering. "
                             "If neither --range nor this is provided, all chromosomes will be plotted.")

    args = parser.parse_args()

    # Call the plotting function
    plot_shared_haplotypes(
        args.haplotype_output_file,
        target_chrom=args.target_chromosome,
        genomic_range=args.range,
        mrca_variant_str=args.mrca_variant,
        gene_name=args.gene_name,
        color_by=args.color_by # Pass the new argument
    )
