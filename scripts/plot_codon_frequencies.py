#!/usr/bin/env python3
"""
Visualize ATG and stop codon frequencies across transcript regions
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[1]

# Load results
df = pd.read_csv(BASE_DIR / 'results' / 'codon_frequencies_by_region.tsv', sep='\t')

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

regions = df['region'].values
x = np.arange(len(regions))
width = 0.35

# Panel 1: Absolute frequencies
ax1.bar(x - width/2, df['atg_per_kb'], width, label='ATG (start)',
        color='steelblue', edgecolor='black', alpha=0.8)
ax1.bar(x + width/2, df['stops_per_kb'], width, label='Stop codons (TAA/TAG/TGA)',
        color='coral', edgecolor='black', alpha=0.8)

ax1.set_ylabel('Codons per kb', fontsize=12, fontweight='bold')
ax1.set_xlabel('Transcript Region', fontsize=12, fontweight='bold')
ax1.set_title('Codon Frequencies by Region', fontsize=13, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(['5\' UTR', 'CDS', '3\' UTR'])
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for i, (atg, stop) in enumerate(zip(df['atg_per_kb'], df['stops_per_kb'])):
    ax1.text(i - width/2, atg + 0.5, f'{atg:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax1.text(i + width/2, stop + 0.5, f'{stop:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

# Panel 2: Fold change relative to CDS
cds_atg = df[df['region'] == 'CDS']['atg_per_kb'].values[0]
cds_stop = df[df['region'] == 'CDS']['stops_per_kb'].values[0]

atg_fold = df['atg_per_kb'] / cds_atg
stop_fold = df['stops_per_kb'] / cds_stop

ax2.bar(x - width/2, atg_fold, width, label='ATG (start)',
        color='steelblue', edgecolor='black', alpha=0.8)
ax2.bar(x + width/2, stop_fold, width, label='Stop codons',
        color='coral', edgecolor='black', alpha=0.8)

ax2.axhline(y=1.0, color='black', linestyle='--', linewidth=2, alpha=0.5, label='CDS baseline')
ax2.set_ylabel('Fold change vs CDS', fontsize=12, fontweight='bold')
ax2.set_xlabel('Transcript Region', fontsize=12, fontweight='bold')
ax2.set_title('Relative Codon Frequencies', fontsize=13, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(['5\' UTR', 'CDS', '3\' UTR'])
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for i, (atg, stop) in enumerate(zip(atg_fold, stop_fold)):
    ax2.text(i - width/2, atg + 0.03, f'{atg:.2f}x', ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax2.text(i + width/2, stop + 0.03, f'{stop:.2f}x', ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.tight_layout()

# Save figure
output_path = BASE_DIR / 'results' / 'figures' / 'codon_frequencies_by_region.png'
Path(output_path).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
print(f"Saved: {output_path}")

plt.show()
