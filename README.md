# SHP — Short Peptide Hunting Pipeline

A Nextflow DSL2 pipeline for discovering short open reading frames (sORFs, 10–100 amino acids) in genomic regions flanking target genes identified by the [allele_search](../allele_search) pipeline.

---

## Biological Motivation

Short ORFs near regulatory genes are a well-documented phenomenon in bacteria. Several classes exist:

- **Translation arrest peptides** (10–30 aa): encoded upstream of protein-localization genes such as *secA*, *secDF*, and *yidC*; when translated, they stall ribosomes and regulate downstream gene expression
- **Leader peptides** (6–20 aa): upstream ORFs in attenuator-based regulation of amino acid biosynthesis operons
- **Toxin-antitoxin sORFs**: often encoded on the strand opposite a toxin or regulatory gene

SHP automates the discovery of these sORFs across large genome collections, starting from the gene-hit coordinates produced by allele_search.

---

## What the Pipeline Does

1. **Extract genomic windows** — for each allele_search gene hit, cuts out a configurable bp window relative to the stop codon (offsets can be negative to overlap the gene body)
2. **Find ORFs** — scans the window for all ORFs ≥ `--min_orf_length` using ATG/GTG/TTG start codons; strand-configurable
3. **HMMER homology search** *(optional)* — translates ORFs to protein and searches against a user-provided query using `phmmer`, `jackhmmer`, or `hmmsearch`; results annotate final outputs
4. **Collect all ORFs** — merges ORFs from all genome × hit combinations into a single FASTA, preserving full genome provenance in headers
5. **Translate and group into alleles** — translates each ORF to protein (GTG/TTG → Met), groups identical proteins into allele families per source gene, and writes output in allele_search allele_info format with per-genome occurrence tables

---

## Requirements

- **Nextflow** (installed on host; not inside Docker)
- **Docker** (recommended) or **conda**
- allele_search output: `blast_summary_filtered/*.csv` (one CSV per target gene)
- A directory of genome FASTA files (one per strain)

---

## Installation

### Docker (recommended)

```bash
docker build -t shp:latest ./docker/
```

### Conda (fallback)

```bash
conda env create -f environment.yml
conda activate shp
```

The conda profile hard-codes an environment path in `nextflow.config`. Update the `beforeScript` line in the `standard` profile if your conda prefix differs.

---

## Quick Start

```bash
nextflow run main.nf \
  --allele_search_dir /path/to/allele_search/results/blast_summary_filtered \
  --genome_dir /path/to/genomes \
  --region_start -50 \
  --region_end 150 \
  --strand_mode both \
  --outdir ./results \
  -profile docker \
  -resume
```

This example searches 200 bp around the stop codon of each target gene hit (50 bp into the gene body, 150 bp downstream) on all reading frames.

---

## HMMER Homology Search (Optional Step 3.5)

Use `--hmmer_query` to search translated ORFs for homologs of a known query. The step is fully optional and backward-compatible — omitting `--hmmer_query` leaves all outputs identical to a run without it.

### Tool selection

| Tool | Query type | Use when |
|---|---|---|
| `phmmer` *(default)* | Single protein FASTA | You have one reference sequence |
| `jackhmmer` | Single protein FASTA | You want iterative refinement (PSI-BLAST-like) |
| `hmmsearch` | HMM profile (`.hmm`) **or** FASTA | You have a curated profile, or multiple reference sequences |

### Examples

```bash
# phmmer: single reference protein


# hmmsearch: pre-built HMM profile
nextflow run main.nf \
  --allele_search_dir ./allele_results/blast_summary_filtered \
  --genome_dir /genomes \
  --region_start 0 --region_end 500 \
  --hmmer_query /data/arrest_peptides.hmm \
  --hmmer_tool hmmsearch \
  --outdir ./results -profile docker -resume

# hmmsearch: build profile from a FASTA of reference sequences (with hmmalign refinement)
nextflow run main.nf \
  --allele_search_dir ./allele_results/blast_summary_filtered \
  --genome_dir /genomes \
  --region_start 0 --region_end 500 \
  --hmmer_query /data/homologs.faa \
  --hmmer_tool hmmsearch \
  --outdir ./results -profile docker -resume
```

### HMM profile building (hmmsearch + FASTA query)

When `--hmmer_tool hmmsearch` is used with a FASTA query (not `.hmm`), the pipeline builds a profile automatically:

1. `hmmbuild` — build an initial HMM from the input FASTA
2. `hmmalign --trim` — re-align the sequences to the initial profile to improve the model
3. `hmmbuild` — rebuild the final HMM from the re-aligned sequences
4. `hmmsearch` — search with the final profile

Set `--hmmer_use_original_alignment true` to skip steps 2–3 and use the input FASTA directly.

---

## Parameters Reference

| Parameter | Default | Description |
|---|---|---|
| `--allele_search_dir` | required | Path to allele_search `blast_summary_filtered/` directory |
| `--genome_dir` | required | Directory of genome FASTA files (`.fasta`, `.fa`, `.fna`, or `.ffn`) |
| `--region_start` | required | bp offset from gene end for window start; negative values overlap the gene body |
| `--region_end` | required | bp offset from gene end for window end; must be > `--region_start` |
| `--strand_mode` | `both` | Reading frames to search: `same_strand`, `other_strand`, or `both` (all 6 frames) |
| `--min_orf_length` | `30` | Minimum ORF length in bp (= 10 aa) |
| `--max_orf_length` | none | Maximum ORF length in bp, excluding stop codon; no upper limit if omitted |
| `--outdir` | `./results` | Output directory |
| `--code_dir` | `${projectDir}/code` | Directory containing Python scripts |
| `--start_codons` | `ATG,GTG,TTG` | Comma-separated start codons to recognize |
| `--batch_size` | `100` | Region FASTAs processed per `FIND_ORFS` task |
| `--extract_batch_size` | `100` | Genomes batched per `EXTRACT_DOWNSTREAM_REGION` task |

### HMMER search parameters (all optional; step is skipped if `--hmmer_query` is not set)

| Parameter | Default | Description |
|---|---|---|
| `--hmmer_query` | none | Path to query file: protein FASTA (`.faa`/`.fasta`) for phmmer/jackhmmer/hmmsearch, or HMM profile (`.hmm`) for hmmsearch |
| `--hmmer_tool` | `phmmer` | HMMER program to use: `phmmer`, `jackhmmer`, or `hmmsearch` |
| `--hmmer_evalue` | `10.0` | Reporting E-value threshold (`-E`) |
| `--hmmer_incE` | `0.001` | Inclusion E-value threshold (`--incE`); controls which hits are considered significant |
| `--hmmer_N` | `5` | Maximum jackhmmer iterations (`-N`); ignored for phmmer/hmmsearch |
| `--hmmer_cpu` | `2` | CPU threads per HMMER process (`--cpu`) |
| `--hmmer_use_original_alignment` | `false` | For hmmsearch + FASTA query: skip the hmmalign re-alignment step and use the input FASTA directly for hmmbuild |

**Strand-aware windowing:** offsets are always relative to the stop codon and mean "downstream past the stop codon" for both forward and reverse strand hits. A `--region_start` of `-50` captures 50 bp inside the gene for both orientations.

---

## Step Modes

Use `--extract` to re-run the full pipeline with updated parameters (e.g. a new window size). Both `--extract` and the default full run execute all pipeline steps.

```bash
nextflow run main.nf --extract \
  --allele_search_dir ./allele_results/blast_summary_filtered \
  --genome_dir /genomes \
  --region_start -50 --region_end 300 \
  --outdir ./results -profile docker -resume
```

---

## Outputs

### `allele_info/{source_gene}/self_reference/`

One subdirectory per target gene from allele_search.

**`allele_info_no_strain_name.csv`** — one row per unique protein allele, sorted by frequency:

| Column | Description |
|---|---|
| `allele_name` | Allele identifier, e.g. `Rgg144_1` (1 = most frequent) |
| `species_count` | Number of unique genomes carrying this allele |
| `species_percentage` | Percentage of total genomes for this source gene |
| `variant_position_count` | `na` |
| `variant_percentage` | `na` |
| `variant_positions(position:reference→changed)` | (empty) |
| `allele_sequence` | Protein sequence |

**`allele_info_with_strain_name.csv`** — same as above, plus `strain_names`: semicolon-separated genome file basenames that carry this allele.

When `--hmmer_query` is set, both CSVs gain two additional trailing columns:

| Column | Description |
|---|---|
| `hmmer_hit_count` | Number of unique genomes where this allele has at least one HMMER hit |
| `hmmer_hit_percentage` | `hmmer_hit_count` / total genomes for this source gene × 100 |

**`allele_sequence_output.fasta`** — one entry per allele:
```
>Rgg144_1|species_count=17|species_percentage=9.44
MSTV...
```

### `per_genome_orfs/{source_gene}.csv`

One row per ORF occurrence (not deduplicated). Links every ORF back to its source genome and assigns its allele name:

| Column | Description |
|---|---|
| `genome` | Genome file basename |
| `contig` | Contig name within that genome |
| `orf_start` | 1-based genomic start (5′ end of ORF) |
| `orf_end` | 1-based genomic end (3′ end including stop codon) |
| `strand` | `+` or `-` |
| `source_gene` | Target gene from allele_search |
| `allele_name` | Allele identifier (links to allele_info CSVs) |
| `protein_sequence` | Translated protein |

When `--hmmer_query` is set, two additional trailing columns are added:

| Column | Description |
|---|---|
| `hmmer_hit` | `True` if this ORF matched the query; `False` otherwise |
| `hmmer_evalue` | Best sequence E-value from HMMER; `NA` if no hit |

### `hmmer/all_hits.tblout` *(only when `--hmmer_query` is set)*

Merged HMMER `--tblout` output across all regions, in standard HMMER tabular format. Columns follow the [HMMER tblout specification](http://hmmer.org/documentation.html): target name, query name, sequence E-value, score, and domain statistics.

### `candidate_orfs/all_orfs_raw.fasta`

Concatenated raw ORF DNA sequences from all genomes and hits, with full header metadata. Intermediate file; input to the allele_info step.

---

## Notes

- **GTG/TTG start codons → Met:** `find_orfs.py` records the actual start codon in the ORF DNA; translation always converts the first codon to Met (bacterial fMet convention).
- **No deduplication at the DNA level:** all ORF occurrences from all genomes are kept, so the same protein can appear in `per_genome_orfs/` multiple times (once per genome that carries it). Deduplication happens at the protein level when building allele groups.
- **Contig boundary clipping:** if the region window extends past a contig boundary, the sequence is silently clipped. Short flanking sequences may yield no ORFs — expected behavior.
- **Docker + symlinked genomes:** symlinks pointing outside the mounted home directory will not resolve inside the container. Use real file paths, or run with `-profile standard` (conda).
- **Work directory size:** the Nextflow `work/` directory grows large over time. Run `nextflow clean -f` periodically.
- **`-resume`:** reuses cached tasks from previous runs. Omit to force a full re-run.

---

## Relationship to allele_search

SHP is a companion to the [allele_search](../allele_search) pipeline. It takes allele_search's `blast_summary_filtered/*.csv` output as its primary input (gene hit coordinates per genome) and writes output in the same allele_info format used by allele_search, making downstream comparison straightforward.
