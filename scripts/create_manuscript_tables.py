#!/usr/bin/env python3
"""Create manuscript-ready supplementary tables from scored overlap tables."""

from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "manuscript_tables"

CDS_CENTRIC_CONTEXTS = {
    "includes_cds",
    "overlaps_upstream",
    "overlaps_downstream",
    "cds_internal",
}

COLUMN_DESCRIPTIONS = [
    ("transcript_id", "GENCODE transcript identifier (ENST ID with version)"),
    ("gene_name", "HGNC gene symbol"),
    ("gene_id", "GENCODE gene identifier (ENSG ID)"),
    ("biotype", "Transcript biotype from GENCODE annotation"),
    ("region_start", "Start position of triple-frame region in transcript coordinates (0-based)"),
    ("region_end", "End position of triple-frame region in transcript coordinates (0-based)"),
    ("region_length_codons", "Length of triple-frame region in codons"),
    ("has_cds_overlap", "Boolean indicating if region overlaps annotated CDS"),
    ("cds_overlap_codons", "Number of codons overlapping with annotated CDS"),
    ("cds_start", "Start position of annotated CDS in transcript coordinates"),
    ("cds_end", "End position of annotated CDS in transcript coordinates"),
    ("frame0_length", "Length of ORF/region in reading frame 0 (codons)"),
    ("frame1_length", "Length of ORF/region in reading frame 1 (codons)"),
    ("frame2_length", "Length of ORF/region in reading frame 2 (codons)"),
    ("rank", "Rank by region length among all triple-frame regions"),
    ("cds_context", "Relationship to annotated CDS"),
    ("total_riboseq_counts", "Total Ribo-seq P-site coverage across the region"),
    ("frame0_riboseq_counts", "Ribo-seq P-site coverage assigned to reading frame 0"),
    ("frame1_riboseq_counts", "Ribo-seq P-site coverage assigned to reading frame 1"),
    ("frame2_riboseq_counts", "Ribo-seq P-site coverage assigned to reading frame 2"),
    ("riboseq_density", "Ribo-seq coverage normalized by region length"),
    ("frame_balance", "Entropy-based balance of coverage across the three frames"),
    ("cds_centric_rank", "Rank by region length among CDS-centric candidates only"),
]


def clean_transcript_id(value):
    """Keep the ENST identifier with version from the FASTA header-derived ID."""
    return str(value).split("|", 1)[0]


def load_table(path):
    df = pd.read_csv(path, sep="\t")
    df["transcript_id"] = df["transcript_id"].map(clean_transcript_id)
    return df


def cds_centric_subset(df):
    subset = df[
        df["cds_context"].isin(CDS_CENTRIC_CONTEXTS)
        & (df["total_riboseq_counts"] > 0)
    ].copy()
    subset = subset.sort_values("region_length_codons", ascending=False).reset_index(drop=True)
    subset["cds_centric_rank"] = range(1, len(subset) + 1)
    return subset


def write_table(df, name):
    out = OUTPUT_DIR / name
    df.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out} ({len(df):,} rows)")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    atg_stop = load_table(RESULTS_DIR / "triple_frame_riboseq" / "atg_stop_with_riboseq.tsv")
    stop_free = load_table(RESULTS_DIR / "triple_frame_riboseq" / "stop_free_with_riboseq.tsv")

    write_table(atg_stop, "unfiltered_atg_stop.tsv")
    write_table(cds_centric_subset(atg_stop), "cds_centric_atg_stop.tsv")
    write_table(cds_centric_subset(stop_free), "cds_centric_stop_free.tsv")

    pd.DataFrame(COLUMN_DESCRIPTIONS, columns=["column_name", "description"]).to_csv(
        OUTPUT_DIR / "COLUMN_DESCRIPTIONS.tsv", sep="\t", index=False
    )
    print(f"Wrote {OUTPUT_DIR / 'COLUMN_DESCRIPTIONS.tsv'}")


if __name__ == "__main__":
    main()
