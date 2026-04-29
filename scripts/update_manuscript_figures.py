#!/usr/bin/env python3
"""
Update manuscript figures based on reviewer feedback:
1. Use transparent bars with only tops colored (not split bars)
2. Change panel labels to lowercase (a, b, c)
3. Increase legend font size
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Rectangle, Patch

# Load data
print("Loading data...")
base_dir = Path(__file__).parent.parent
results_dir = base_dir / 'results' / 'triple_frame_riboseq'
revcomp_dir = base_dir / 'results' / 'revcomp_triple_frame'

df_stop_free = pd.read_csv(results_dir / 'stop_free_with_riboseq.tsv', sep='\t')
df_atg_stop = pd.read_csv(results_dir / 'atg_stop_with_riboseq.tsv', sep='\t')
df_rc_stop_free = pd.read_csv(revcomp_dir / 'revcomp_stop_free_triple_frame.tsv', sep='\t')
df_rc_atg_stop = pd.read_csv(revcomp_dir / 'revcomp_atg_stop_triple_frame.tsv', sep='\t')

# Filter CDS-centric
def filter_cds_centric(df):
    cds_centric_contexts = ['includes_cds', 'overlaps_upstream', 'overlaps_downstream', 'cds_internal']
    filtered = df[
        (df['cds_context'].isin(cds_centric_contexts)) &
        (df['total_riboseq_counts'] > 0)
    ].copy()
    return filtered

cds_atg_stop = filter_cds_centric(df_atg_stop)

print(f"Loaded {len(cds_atg_stop):,} CDS-centric ATG-to-STOP regions")
print(f"Loaded {len(df_stop_free):,} unfiltered Stop-Free regions")
print(f"Loaded {len(df_atg_stop):,} unfiltered ATG-to-STOP regions")


def create_main_figure_updated(df, df_revcomp, method_name, filename_suffix, color):
    """Create 1x3 main figure with transparent bars and colored tops"""

    if len(df) == 0:
        print(f"No candidates found for {method_name}!")
        return

    # Create figure
    fig = plt.figure(figsize=(21, 6))
    gs = fig.add_gridspec(1, 3, wspace=0.3)

    # Bins
    bin_width = 10
    bins = np.arange(df['region_length_codons'].min(), df['region_length_codons'].max() + bin_width, bin_width)
    max_len_c = max(df['region_length_codons'].max(), df_revcomp['region_length_codons'].max())
    bins_c = np.arange(30, max_len_c + 10, 10)

    # Get SRD5A1
    srd5a1 = df[df['gene_name'] == 'SRD5A1']
    srd5a1_length = srd5a1.iloc[0]['region_length_codons'] if len(srd5a1) > 0 else None
    srd5a1_row = srd5a1.iloc[0] if len(srd5a1) > 0 else None

    # ==================================================================
    # PANEL a: Ribo-Seq Support (transparent bars with colored tops)
    # ==================================================================
    ax1 = fig.add_subplot(gs[0, 0])

    # Create Ribo-Seq tiers
    df_copy = df.copy()
    df_copy['riboseq_tier'] = pd.cut(df_copy['riboseq_density'],
                                      bins=[0, 1, 100, 1000, np.inf],
                                      labels=['None', 'Low (1-100)', 'Medium (100-1K)', 'High (>1K)'])

    tier_colors = {
        'None': '#d3d3d3',
        'Low (1-100)': '#fee5d9',
        'Medium (100-1K)': '#fc9272',
        'High (>1K)': '#de2d26'
    }

    # Calculate stacked heights to know where to draw colored tops
    tier_heights = {}
    bottom = np.zeros(len(bins)-1)

    for tier in ['None', 'Low (1-100)', 'Medium (100-1K)', 'High (>1K)']:
        tier_data = df_copy[df_copy['riboseq_tier'] == tier]['region_length_codons']
        if len(tier_data) > 0:
            counts, _ = np.histogram(tier_data, bins=bins)
            tier_heights[tier] = (bottom.copy(), counts.copy())
            bottom += counts

    # Draw transparent bars for full height
    total_counts, _ = np.histogram(df['region_length_codons'], bins=bins)
    ax1.bar(bins[:-1], total_counts, width=bin_width,
           color='none', edgecolor='black', linewidth=0.5)

    # Draw colored tops for each tier
    top_height = 0.05  # Height of colored top in data coordinates (will scale with log)

    for tier in ['None', 'Low (1-100)', 'Medium (100-1K)', 'High (>1K)']:
        if tier in tier_heights:
            tier_bottom, tier_counts = tier_heights[tier]
            top_positions = tier_bottom + tier_counts

            # For each bin, draw a colored rectangle at the top
            for i, (x, top_y) in enumerate(zip(bins[:-1], top_positions)):
                if tier_counts[i] > 0:
                    # Calculate height for colored top (percentage of bar height)
                    # Use a fixed percentage for log scale
                    colored_height = max(top_y * 0.15, 1)  # 15% of bar height or minimum 1

                    rect = Rectangle((x, top_y - colored_height), bin_width, colored_height,
                                   facecolor=tier_colors[tier], edgecolor='none',
                                   alpha=0.9, label=tier if i == 0 and tier_counts[i] > 0 else "")
                    ax1.add_patch(rect)

    # Manually create legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=tier_colors[tier], label=tier, alpha=0.9)
                      for tier in ['None', 'Low (1-100)', 'Medium (100-1K)', 'High (>1K)']]

    if srd5a1_length:
        ax1.axvline(srd5a1_length, color='blue', linestyle='--', linewidth=2.5,
                   label=f'SRD5A1 ({srd5a1_row["riboseq_density"]:.0f} rpf/codon)')
        legend_elements.append(plt.Line2D([0], [0], color='blue', linestyle='--', linewidth=2.5,
                                         label=f'SRD5A1 ({srd5a1_row["riboseq_density"]:.0f} rpf/codon)'))

    ax1.set_xlabel('Region length (codons)', fontsize=12)
    ax1.set_ylabel('Number of regions', fontsize=12)
    ax1.set_yscale('log')
    ax1.set_title('a. Ribo-Seq Support', fontsize=14, fontweight='bold', loc='left')
    ax1.legend(handles=legend_elements, fontsize=10, loc='upper right')
    ax1.grid(True, alpha=0.3, axis='y')

    # ==================================================================
    # PANEL b: CDS Context (transparent bars with colored tops)
    # ==================================================================
    ax2 = fig.add_subplot(gs[0, 1])

    context_colors = {
        'includes_cds': '#2ecc71',
        'overlaps_upstream': '#3498db',
        'overlaps_downstream': '#9b59b6',
        'cds_internal': '#e74c3c',
        'no_overlap': '#bdc3c7',
        'no_cds_annotation': '#7f8c8d'
    }

    # Calculate stacked heights
    context_heights = {}
    bottom = np.zeros(len(bins)-1)

    for context in ['includes_cds', 'overlaps_upstream', 'overlaps_downstream', 'cds_internal', 'no_overlap', 'no_cds_annotation']:
        context_data = df[df['cds_context'] == context]['region_length_codons']
        if len(context_data) > 0:
            counts, _ = np.histogram(context_data, bins=bins)
            context_heights[context] = (bottom.copy(), counts.copy())
            bottom += counts

    # Draw transparent bars
    total_counts, _ = np.histogram(df['region_length_codons'], bins=bins)
    ax2.bar(bins[:-1], total_counts, width=bin_width,
           color='none', edgecolor='black', linewidth=0.5)

    # Draw colored tops
    for context in ['includes_cds', 'overlaps_upstream', 'overlaps_downstream', 'cds_internal', 'no_overlap', 'no_cds_annotation']:
        if context in context_heights:
            context_bottom, context_counts = context_heights[context]
            top_positions = context_bottom + context_counts

            for i, (x, top_y) in enumerate(zip(bins[:-1], top_positions)):
                if context_counts[i] > 0:
                    colored_height = max(top_y * 0.15, 1)

                    rect = Rectangle((x, top_y - colored_height), bin_width, colored_height,
                                   facecolor=context_colors[context], edgecolor='none',
                                   alpha=0.9)
                    ax2.add_patch(rect)

    # Create legend
    legend_elements = [Patch(facecolor=context_colors[ctx], label=ctx.replace('_', ' '), alpha=0.9)
                      for ctx in ['includes_cds', 'overlaps_upstream', 'overlaps_downstream', 'cds_internal', 'no_overlap', 'no_cds_annotation']
                      if ctx in context_heights]

    if srd5a1_length:
        ax2.axvline(srd5a1_length, color='blue', linestyle='--', linewidth=2.5)
        legend_elements.append(plt.Line2D([0], [0], color='blue', linestyle='--', linewidth=2.5,
                                         label=f'SRD5A1 ({srd5a1_row["cds_context"]})'))

    ax2.set_xlabel('Region length (codons)', fontsize=12)
    ax2.set_ylabel('Number of regions', fontsize=12)
    ax2.set_yscale('log')
    ax2.set_title('b. CDS Context', fontsize=14, fontweight='bold', loc='left')
    ax2.legend(handles=legend_elements, fontsize=10, loc='upper right')
    ax2.grid(True, alpha=0.3, axis='y')

    # ==================================================================
    # PANEL c: Control (overlapping histograms - use same transparent approach)
    # ==================================================================
    ax3 = fig.add_subplot(gs[0, 2])

    fwd_lengths = df['region_length_codons']
    rc_lengths = df_revcomp['region_length_codons']

    # Get histogram values
    fwd_counts, _ = np.histogram(fwd_lengths, bins=bins_c)
    rc_counts, _ = np.histogram(rc_lengths, bins=bins_c)

    # Draw transparent bars with colored tops
    # Forward strand
    ax3.bar(bins_c[:-1], fwd_counts, width=10,
           color='none', edgecolor='black', linewidth=0.5, alpha=0.8)

    for i, (x, y) in enumerate(zip(bins_c[:-1], fwd_counts)):
        if y > 0:
            colored_height = max(y * 0.15, 1)
            rect = Rectangle((x, y - colored_height), 10, colored_height,
                           facecolor=color, edgecolor='none', alpha=0.7)
            ax3.add_patch(rect)

    # Reverse complement (draw on same bins, will overlap)
    ax3.bar(bins_c[:-1], rc_counts, width=10,
           color='none', edgecolor='black', linewidth=0.5, alpha=0.8)

    for i, (x, y) in enumerate(zip(bins_c[:-1], rc_counts)):
        if y > 0:
            colored_height = max(y * 0.15, 1)
            rect = Rectangle((x, y - colored_height), 10, colored_height,
                           facecolor='gray', edgecolor='none', alpha=0.7)
            ax3.add_patch(rect)

    # Create legend
    legend_elements = [
        Patch(facecolor=color, label=f'Forward (n={len(fwd_lengths):,})', alpha=0.7),
        Patch(facecolor='gray', label=f'RevComp (n={len(rc_lengths):,})', alpha=0.7)
    ]

    ax3.set_xlabel('Region length (codons)', fontsize=12)
    ax3.set_ylabel('Frequency', fontsize=12)
    ax3.set_yscale('log')
    ax3.set_title('c. Control Comparison', fontsize=14, fontweight='bold', loc='left')
    ax3.legend(handles=legend_elements, fontsize=10, loc='upper right')
    ax3.grid(True, alpha=0.3, axis='y')

    # Save
    fig_path = base_dir / 'results' / 'figures' / f'manuscript_{filename_suffix}.png'
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(fig_path).replace('.png', '.pdf'), bbox_inches='tight')
    print(f"Saved: {fig_path}")
    plt.close()


def create_supplementary_figure_updated(df_stop_free, df_rc_stop_free, df_atg_stop, df_rc_atg_stop):
    """Create 2x3 supplementary with lowercase labels and transparent bars"""

    fig = plt.figure(figsize=(21, 12))
    gs = fig.add_gridspec(2, 3, wspace=0.3, hspace=0.3)

    datasets = [
        (df_stop_free, df_rc_stop_free, 'Stop-Free', 'steelblue', 0),
        (df_atg_stop, df_rc_atg_stop, 'ATG-to-STOP', 'forestgreen', 1)
    ]

    panel_labels = ['a', 'b', 'c', 'd', 'e', 'f']

    for row_idx, (df, df_revcomp, method_name, color, _) in enumerate(datasets):
        bin_width = 10
        bins = np.arange(df['region_length_codons'].min(), df['region_length_codons'].max() + bin_width, bin_width)
        max_len_c = max(df['region_length_codons'].max(), df_revcomp['region_length_codons'].max())
        bins_c = np.arange(30, max_len_c + 10, 10)

        srd5a1 = df[df['gene_name'] == 'SRD5A1']
        srd5a1_length = srd5a1.iloc[0]['region_length_codons'] if len(srd5a1) > 0 else None
        srd5a1_row = srd5a1.iloc[0] if len(srd5a1) > 0 else None

        # PANEL 1: Ribo-Seq
        ax1 = fig.add_subplot(gs[row_idx, 0])

        df_copy = df.copy()
        df_copy['riboseq_tier'] = pd.cut(df_copy['riboseq_density'],
                                          bins=[0, 1, 100, 1000, np.inf],
                                          labels=['None', 'Low (1-100)', 'Medium (100-1K)', 'High (>1K)'])

        tier_colors = {'None': '#d3d3d3', 'Low (1-100)': '#fee5d9', 'Medium (100-1K)': '#fc9272', 'High (>1K)': '#de2d26'}

        # Stack heights
        tier_heights = {}
        bottom = np.zeros(len(bins)-1)
        for tier in ['None', 'Low (1-100)', 'Medium (100-1K)', 'High (>1K)']:
            tier_data = df_copy[df_copy['riboseq_tier'] == tier]['region_length_codons']
            if len(tier_data) > 0:
                counts, _ = np.histogram(tier_data, bins=bins)
                tier_heights[tier] = (bottom.copy(), counts.copy())
                bottom += counts

        # Transparent bars
        total_counts, _ = np.histogram(df['region_length_codons'], bins=bins)
        ax1.bar(bins[:-1], total_counts, width=bin_width, color='none', edgecolor='black', linewidth=0.5)

        # Colored tops
        for tier in ['None', 'Low (1-100)', 'Medium (100-1K)', 'High (>1K)']:
            if tier in tier_heights:
                tier_bottom, tier_counts = tier_heights[tier]
                top_positions = tier_bottom + tier_counts
                for i, (x, top_y) in enumerate(zip(bins[:-1], top_positions)):
                    if tier_counts[i] > 0:
                        colored_height = max(top_y * 0.15, 1)
                        rect = Rectangle((x, top_y - colored_height), bin_width, colored_height,
                                       facecolor=tier_colors[tier], edgecolor='none', alpha=0.9)
                        ax1.add_patch(rect)

        legend_elements = [Patch(facecolor=tier_colors[tier], label=tier, alpha=0.9)
                          for tier in ['None', 'Low (1-100)', 'Medium (100-1K)', 'High (>1K)'] if tier in tier_heights]
        if srd5a1_length:
            ax1.axvline(srd5a1_length, color='blue', linestyle='--', linewidth=2.5)
            legend_elements.append(plt.Line2D([0], [0], color='blue', linestyle='--', linewidth=2.5, label='SRD5A1'))

        ax1.set_xlabel('Region length (codons)', fontsize=11)
        ax1.set_ylabel('Number of regions', fontsize=11)
        ax1.set_yscale('log')
        ax1.set_title(f'{panel_labels[row_idx*3]}. {method_name} - Ribo-Seq', fontsize=13, fontweight='bold', loc='left')
        ax1.legend(handles=legend_elements, fontsize=10, loc='upper right')
        ax1.grid(True, alpha=0.3, axis='y')

        # PANEL 2: CDS Context
        ax2 = fig.add_subplot(gs[row_idx, 1])

        context_colors = {
            'includes_cds': '#2ecc71', 'overlaps_upstream': '#3498db', 'overlaps_downstream': '#9b59b6',
            'cds_internal': '#e74c3c', 'no_overlap': '#bdc3c7', 'no_cds_annotation': '#7f8c8d'
        }

        context_heights = {}
        bottom = np.zeros(len(bins)-1)
        for context in ['includes_cds', 'overlaps_upstream', 'overlaps_downstream', 'cds_internal', 'no_overlap', 'no_cds_annotation']:
            context_data = df[df['cds_context'] == context]['region_length_codons']
            if len(context_data) > 0:
                counts, _ = np.histogram(context_data, bins=bins)
                context_heights[context] = (bottom.copy(), counts.copy())
                bottom += counts

        ax2.bar(bins[:-1], total_counts, width=bin_width, color='none', edgecolor='black', linewidth=0.5)

        for context in ['includes_cds', 'overlaps_upstream', 'overlaps_downstream', 'cds_internal', 'no_overlap', 'no_cds_annotation']:
            if context in context_heights:
                context_bottom, context_counts = context_heights[context]
                top_positions = context_bottom + context_counts
                for i, (x, top_y) in enumerate(zip(bins[:-1], top_positions)):
                    if context_counts[i] > 0:
                        colored_height = max(top_y * 0.15, 1)
                        rect = Rectangle((x, top_y - colored_height), bin_width, colored_height,
                                       facecolor=context_colors[context], edgecolor='none', alpha=0.9)
                        ax2.add_patch(rect)

        legend_elements = [Patch(facecolor=context_colors[ctx], label=ctx.replace('_', ' '), alpha=0.9)
                          for ctx in ['includes_cds', 'overlaps_upstream', 'overlaps_downstream', 'cds_internal', 'no_overlap', 'no_cds_annotation']
                          if ctx in context_heights]
        if srd5a1_length:
            ax2.axvline(srd5a1_length, color='blue', linestyle='--', linewidth=2.5)
            legend_elements.append(plt.Line2D([0], [0], color='blue', linestyle='--', linewidth=2.5, label='SRD5A1'))

        ax2.set_xlabel('Region length (codons)', fontsize=11)
        ax2.set_ylabel('Number of regions', fontsize=11)
        ax2.set_yscale('log')
        ax2.set_title(f'{panel_labels[row_idx*3+1]}. {method_name} - CDS Context', fontsize=13, fontweight='bold', loc='left')
        ax2.legend(handles=legend_elements, fontsize=10, loc='upper right')
        ax2.grid(True, alpha=0.3, axis='y')

        # PANEL 3: Control
        ax3 = fig.add_subplot(gs[row_idx, 2])

        fwd_lengths = df['region_length_codons']
        rc_lengths = df_revcomp['region_length_codons']

        fwd_counts, _ = np.histogram(fwd_lengths, bins=bins_c)
        rc_counts, _ = np.histogram(rc_lengths, bins=bins_c)

        ax3.bar(bins_c[:-1], fwd_counts, width=10, color='none', edgecolor='black', linewidth=0.5)
        for i, (x, y) in enumerate(zip(bins_c[:-1], fwd_counts)):
            if y > 0:
                colored_height = max(y * 0.15, 1)
                rect = Rectangle((x, y - colored_height), 10, colored_height,
                               facecolor=color, edgecolor='none', alpha=0.7)
                ax3.add_patch(rect)

        ax3.bar(bins_c[:-1], rc_counts, width=10, color='none', edgecolor='black', linewidth=0.5)
        for i, (x, y) in enumerate(zip(bins_c[:-1], rc_counts)):
            if y > 0:
                colored_height = max(y * 0.15, 1)
                rect = Rectangle((x, y - colored_height), 10, colored_height,
                               facecolor='gray', edgecolor='none', alpha=0.7)
                ax3.add_patch(rect)

        legend_elements = [
            Patch(facecolor=color, label=f'Forward (n={len(fwd_lengths):,})', alpha=0.7),
            Patch(facecolor='gray', label=f'RevComp (n={len(rc_lengths):,})', alpha=0.7)
        ]

        ax3.set_xlabel('Region length (codons)', fontsize=11)
        ax3.set_ylabel('Frequency', fontsize=11)
        ax3.set_yscale('log')
        ax3.set_title(f'{panel_labels[row_idx*3+2]}. {method_name} - Control', fontsize=13, fontweight='bold', loc='left')
        ax3.legend(handles=legend_elements, fontsize=10, loc='upper right')
        ax3.grid(True, alpha=0.3, axis='y')

    fig_path = base_dir / 'results' / 'figures' / 'manuscript_supplementary.png'
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.savefig(str(fig_path).replace('.png', '.pdf'), bbox_inches='tight')
    print(f"Saved: {fig_path}")
    plt.close()


# Generate figures
print("\n" + "="*70)
print("GENERATING UPDATED FIGURES")
print("="*70)

print("\nCreating main figure (ATG-to-STOP, CDS-filtered)...")
create_main_figure_updated(cds_atg_stop, df_rc_atg_stop, 'ATG-to-STOP (Main)', 'atg_stop', 'forestgreen')

print("\nCreating supplementary figure (2x3)...")
create_supplementary_figure_updated(df_stop_free, df_rc_stop_free, df_atg_stop, df_rc_atg_stop)

print("\n" + "="*70)
print("DONE! Updated figures saved with:")
print("  ✓ Transparent bars with colored tops")
print("  ✓ Lowercase panel labels (a, b, c)")
print("  ✓ Increased legend font size (10pt)")
print("="*70)
