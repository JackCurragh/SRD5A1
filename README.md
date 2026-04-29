# SRD5A1 Triple-Frame Translation Analysis

Source code for the genome-wide analyses reported in:

Human SRD5A1 as a case of gene expression indel-resistance in triple-coding region.

This repository contains the notebooks and helper scripts used to identify regions of human transcripts with overlapping open reading frames in all three reading frames, classify their CDS context, score ribosome profiling support, and generate manuscript figures/tables.

## Repository Contents

- `notebooks/01_triple_frame_regions.ipynb`: scans GENCODE v45 transcripts for triple-frame regions using stop-free and ATG-to-stop definitions.
- `notebooks/02_triple_frame_cds_context.ipynb`: classifies candidate regions by position relative to annotated CDSs.
- `notebooks/03_riboseq_scoring.ipynb`: maps transcript-space candidates to genomic coordinates and scores Ribo-seq bigWig coverage by frame.
- `scripts/create_manuscript_tables.py`: creates manuscript-ready supplementary tables from the scored overlap tables.
- `scripts/update_manuscript_figures.py`: regenerates the manuscript distribution figures from processed tables.
- `scripts/analyze_codon_frequencies.py` and `scripts/plot_codon_frequencies.py`: supporting codon-frequency analysis.
- `pixi.toml` and `pixi.lock`: reproducible software environment.

Large input files and derived genome-wide tables are not version-controlled in this source repository. They are generated from public GENCODE v45 annotations, transcript sequences, and RiboCrypt/RiboSeq.org bigWig resources described in the manuscript.

## Required Inputs

Place the following files in the repository root before running the notebooks:

- `annotation.gtf`: GENCODE v45 human GTF annotation.
- `transcripts.fa`: GENCODE v45 human transcript FASTA.
- `TransCODE_GENCODE_v45_RiboSeqORFs_RiboCrypt_all_forward_unique.bigWig`: forward-strand RiboCrypt/RiboSeq.org P-site coverage.
- `TransCODE_GENCODE_v45_RiboSeqORFs_RiboCrypt_all_reverse_unique.bigWig`: reverse-strand RiboCrypt/RiboSeq.org P-site coverage.

The RiboCrypt bigWig files are available from the public EBI FTP location cited in the manuscript.

## Environment

Install dependencies with Pixi:

```bash
pixi install
```

Then start Jupyter from the Pixi environment:

```bash
pixi run jupyter notebook
```

The main dependencies are Python, pandas, pybedtools, intervaltree, pyBigWig, Biopython, matplotlib, and seaborn.

## Workflow

Run the notebooks in order:

1. `notebooks/01_triple_frame_regions.ipynb`
2. `notebooks/02_triple_frame_cds_context.ipynb`
3. `notebooks/03_riboseq_scoring.ipynb`

The figure script expects processed tables under `results/triple_frame_riboseq/` and `results/revcomp_triple_frame/`:

```bash
python scripts/create_manuscript_tables.py
python scripts/update_manuscript_figures.py
```

Supporting codon-frequency analysis:

```bash
python scripts/analyze_codon_frequencies.py
python scripts/plot_codon_frequencies.py
```

## License

This source code is released under the MIT License. See `LICENSE`.

## Citation

Please cite the archived software DOI for the exact version used in the manuscript. Citation metadata for this repository is provided in `CITATION.cff`.
