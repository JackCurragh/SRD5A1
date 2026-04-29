#!/usr/bin/env python3
"""
Analyze ATG and stop codon frequencies in different transcript regions
to understand why 3' UTRs have more triple-frame regions.
"""

import pandas as pd
from pathlib import Path
from collections import defaultdict

BASE_DIR = Path(__file__).resolve().parents[1]

print("Loading data...")
df = pd.read_csv(BASE_DIR / 'results' / 'manuscript_tables' / 'unfiltered_atg_stop.tsv', sep='\t')
print(f"Loaded {len(df):,} regions")

# Load sequences
print("\nLoading transcript sequences...")
sequences = {}
with open(BASE_DIR / 'transcripts.fa') as f:
    current_id = None
    current_seq = []
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_id:
                sequences[current_id] = ''.join(current_seq)
            current_id = line[1:].split('|')[0]  # Get ENST ID
            current_seq = []
        else:
            current_seq.append(line.upper())
    if current_id:
        sequences[current_id] = ''.join(current_seq)

print(f"Loaded {len(sequences):,} transcript sequences")

# Initialize statistics
stats = defaultdict(lambda: {'total_nt': 0, 'atg': 0, 'taa': 0, 'tag': 0, 'tga': 0})

print("\nAnalyzing codon frequencies by region type...")

# Get transcripts with CDS annotations and triple-frame regions
transcripts_to_analyze = df[
    (df['cds_context'] != 'no_cds_annotation') &
    (df['cds_start'].notna())
]['transcript_id'].unique()

print(f"Analyzing {len(transcripts_to_analyze):,} unique transcripts with CDS annotations")

analyzed = 0
for transcript_id_full in transcripts_to_analyze:
    # Extract just the ENST ID from the full transcript name
    transcript_id = transcript_id_full.split('|')[0]

    if transcript_id not in sequences:
        continue

    seq = sequences[transcript_id]

    # Get CDS coordinates from any region in this transcript
    transcript_regions = df[df['transcript_id'] == transcript_id_full]
    if len(transcript_regions) == 0:
        continue

    cds_start = int(transcript_regions.iloc[0]['cds_start'])
    cds_end = int(transcript_regions.iloc[0]['cds_end'])

    # Extract regions (using 0-based indexing)
    utr5 = seq[:cds_start]
    cds = seq[cds_start:cds_end]
    utr3 = seq[cds_end:]

    # Count codons in each region
    for region_name, region_seq in [('5_UTR', utr5), ('CDS', cds), ('3_UTR', utr3)]:
        if len(region_seq) == 0:
            continue

        stats[region_name]['total_nt'] += len(region_seq)
        stats[region_name]['atg'] += region_seq.count('ATG')
        stats[region_name]['taa'] += region_seq.count('TAA')
        stats[region_name]['tag'] += region_seq.count('TAG')
        stats[region_name]['tga'] += region_seq.count('TGA')

    analyzed += 1
    if analyzed % 1000 == 0:
        print(f"  Processed {analyzed:,} transcripts...", end='\r')

print(f"\n\nAnalyzed {analyzed:,} transcripts")

# Calculate frequencies (per kb)
print("\n" + "="*80)
print("CODON FREQUENCIES (per kb of sequence)")
print("="*80)

results = []
for region in ['5_UTR', 'CDS', '3_UTR']:
    s = stats[region]
    total_kb = s['total_nt'] / 1000

    if total_kb == 0:
        continue

    atg_per_kb = s['atg'] / total_kb
    total_stops = s['taa'] + s['tag'] + s['tga']
    stops_per_kb = total_stops / total_kb

    results.append({
        'region': region,
        'total_nt': s['total_nt'],
        'atg_count': s['atg'],
        'stop_count': total_stops,
        'atg_per_kb': atg_per_kb,
        'stops_per_kb': stops_per_kb,
        'atg_stop_ratio': atg_per_kb / stops_per_kb if stops_per_kb > 0 else 0
    })

    print(f"\n{region}:")
    print(f"  Total nucleotides: {s['total_nt']:,} ({total_kb:.1f} kb)")
    print(f"  ATG count: {s['atg']:,} ({atg_per_kb:.2f} per kb)")
    print(f"  Stop codons: {total_stops:,} ({stops_per_kb:.2f} per kb)")
    print(f"    TAA: {s['taa']:,}")
    print(f"    TAG: {s['tag']:,}")
    print(f"    TGA: {s['tga']:,}")
    print(f"  ATG/Stop ratio: {atg_per_kb / stops_per_kb:.3f}")

# Summary comparison
print("\n" + "="*80)
print("RELATIVE TO CDS")
print("="*80)

cds_result = next(r for r in results if r['region'] == 'CDS')

for result in results:
    if result['region'] == 'CDS':
        continue

    atg_fold = result['atg_per_kb'] / cds_result['atg_per_kb']
    stop_fold = result['stops_per_kb'] / cds_result['stops_per_kb']

    print(f"\n{result['region']} vs CDS:")
    print(f"  ATG frequency: {atg_fold:.2f}x")
    print(f"  Stop frequency: {stop_fold:.2f}x")
    print(f"  → {'MORE' if atg_fold > 1 else 'FEWER'} ATGs ({abs(1-atg_fold)*100:.1f}% {'higher' if atg_fold > 1 else 'lower'})")
    print(f"  → {'MORE' if stop_fold > 1 else 'FEWER'} stops ({abs(1-stop_fold)*100:.1f}% {'higher' if stop_fold > 1 else 'lower'})")

# Save results
df_results = pd.DataFrame(results)
output_path = BASE_DIR / 'results' / 'codon_frequencies_by_region.tsv'
df_results.to_csv(output_path, sep='\t', index=False)
print(f"\n\nResults saved to: {output_path}")
