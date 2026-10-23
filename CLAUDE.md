# SHP (Short Peptide Hunting) Pipeline — CLAUDE.md

## Project Context

This is a bioinformatics project using Python. When debugging, always check input data constraints (e.g., sequence lengths, file formats) before diving into code logic.

## Workflow

When working through a checklist or TODO list, always re-read the checklist file after completing each major item to track progress and avoid missing steps.

---

## Project Purpose

Nextflow DSL2 pipeline for discovering short ORFs (sORFs, 10–100 aa) in the genomic region flanking target genes found by the **allele_search** pipeline (`../allele_search`). Given allele_search's `blast_summary_filtered/*.csv` output and a directory of genome FASTA files, SHP:

1. Extracts a configurable genomic window relative to each target gene hit (by bp offset from gene end; window bounds can be negative, overlapping with the gene itself)
2. Finds all ORFs ≥ `min_orf_length` within that window (on same strand, opposite strand, or both)
3. Deduplicates candidate ORFs and translates them to protein (synonymous DNA variants collapse to one protein)
4. Clusters the translated proteins with CD-HIT into protein families

**Output:** A representative protein FASTA (`all_proteins_cdhit.fasta`) and cluster membership file (`all_proteins_cdhit.fasta.clstr`) for each run. Use `code/split_clusters.py` to split clusters into per-cluster FASTAs for downstream analysis.

**Biological motivation:** Short ORFs near regulatory genes are a validated phenomenon in bacteria:
- **Translation arrest peptides** (10–30 aa): encoded upstream of protein-localization genes (*secA*, *secDF*, *yidC*); stall ribosomes to regulate downstream expression
- **Leader peptides** (6–20 aa): upstream ORFs in attenuator-based regulation (e.g., amino acid biosynthesis operons)
- **Toxin-antitoxin sORFs**: often encoded on the strand opposite of a toxin or regulatory gene

The "downstream" window is defined spatially (by genomic coordinates), not biologically — it can be on either strand of the allele_search hit.

**Relationship to allele_search:** SHP uses allele_search output as input. The core logic is ORF finding and protein clustering; downstream allele analysis is done separately.

---

## Tech Stack

| Tool | Version | Role |
|---|---|---|
| Nextflow | DSL2 | Workflow orchestration |
| Python | 3.9 | All custom scripts (`code/`) |
| BioPython | 1.85 | Sequence I/O |
| CD-HIT | latest | Protein clustering |
| Docker | latest | Dependency containerization |

**Docker image:** `shp:latest` (defined in `docker/Dockerfile`)
**Conda fallback:** `shp` environment

---

## Environment Setup

### Docker (recommended)
```bash
# Build the image once
docker build -t shp:latest ./docker/

# Run with Docker profile
nextflow run main.nf --help -profile docker
```

### Conda (fallback)
```bash
conda env create -f environment.yml
conda activate shp
nextflow run main.nf --help
```

`nextflow.config` hard-codes the conda env path in the `standard` profile's `beforeScript`. **Update this path if the conda prefix changes.**

---

## Running the Pipeline

### Full pipeline (most common)
```bash
nextflow run main.nf \
  --allele_search_dir /path/to/allele_search/results/blast_summary_filtered \
  --genome_dir /path/to/genomes \
  --region_start 0 \
  --region_end 500 \
  --strand_mode both \
  --outdir ./results \
  -profile docker \
  -resume
```

### Step-specific modes
Each mode reads from and writes to `--outdir`. Prior steps must have completed first.

```bash
# Re-run region extraction + ORF finding + translation only
nextflow run main.nf --extract \
  --allele_search_dir ./allele_results/blast_summary_filtered \
  --genome_dir /genomes --region_start -50 --region_end 300 \
  --outdir ./results -profile docker

# Re-run CD-HIT clustering only (reuses outdir/candidate_orfs/)
nextflow run main.nf --cluster \
  --outdir ./results -profile docker
```

Alternative syntax: `--step extract|cluster`

### Splitting clusters after the pipeline
```bash
python code/split_clusters.py \
  results/cdhit/all_proteins_cdhit.fasta.clstr \
  results/candidate_orfs/candidate_orfs_protein.fasta \
  -o results/clusters/
```

---

## Directory & Module Map

```
main.nf                          # Workflow entrypoint; routes step modes
nextflow.config                  # All parameters with defaults; Docker + conda profiles
docker/
  Dockerfile                     # Single image: Python 3.9, CD-HIT
modules/
  extract_region.nf              # PARSE_ALLELE_INPUT, EXTRACT_DOWNSTREAM_REGION
  orf_find.nf                    # FIND_ORFS, COLLECT_ORFS
  translate.nf                   # TRANSLATE_CANDIDATES (DNA → protein, deduplicate)
  cluster.nf                     # CD_HIT_CLUSTER
code/
  header_utils.py                # Canonical FASTA header parser/builder (shared utility)
  extract_downstream.py          # Extract genome region by coord + strand
  find_orfs.py                   # ORF finder: configurable strand, min length
  translate_candidates.py        # Translate DNA ORFs to protein; dedup synonymous variants
  split_clusters.py              # Standalone: split CD-HIT clusters into per-cluster FASTAs
results/                         # Pipeline outputs (gitignored)
test/                            # Test data and documentation
```

---

## Pipeline Steps

| Step | Process | Module | Description |
|---|---|---|---|
| 1 | `PARSE_ALLELE_INPUT` | extract_region.nf | Read allele_search `blast_summary_filtered/*.csv` → per-genome gene coordinates (contig, SubjectStart, SubjectEnd) |
| 2 | `EXTRACT_DOWNSTREAM_REGION` | extract_region.nf | For each genome × hit, extract `[gene_end + region_start, gene_end + region_end]` bp window from genome FASTA; strand-aware; clips to contig boundaries |
| 3 | `FIND_ORFS` | orf_find.nf | Find all ORFs ≥ `min_orf_length` in each extracted region; reads frames according to `strand_mode`; ATG+GTG+TTG start codons |
| 4 | `COLLECT_ORFS` | orf_find.nf | Merge ORFs from all hits; deduplicate by exact DNA sequence → candidate ORF FASTA |
| 5 | `TRANSLATE_CANDIDATES` | translate.nf | Translate DNA ORFs to protein in frame 0; GTG/TTG start → Met; deduplicate synonymous variants |
| 6 | `CD_HIT_CLUSTER` | cluster.nf | Cluster proteins by sequence identity; output representative FASTA + .clstr membership file |

---

## Key Conventions & Design Decisions

### FASTA header format (internal)
All scripts communicate metadata through structured FASTA headers, parsed and built exclusively via `header_utils.py`. **Do not change this format without updating `header_utils.py`.**
```
>ID|contig=NAME|E-value=VALUE|AlignLen=LEN|SeqIdentity=%|SubjectStart=POS|SubjectEnd=POS|source_gene=GENE
```
The `source_gene` field (new in SHP; absent in allele_search) tracks which allele_search target gene produced each ORF candidate.

### Strand detection (same convention as allele_search)
**SubjectStart > SubjectEnd signals reverse strand.** Used in `extract_downstream.py` and `find_orfs.py`. No explicit strand column — strand is inferred from coordinates.

### Strand-aware gene end calculation
`gene_end` is always the **stop codon position** (3′ end of the transcript):
- **Forward strand** (SubjectStart ≤ SubjectEnd): gene end = SubjectEnd (higher coord)
- **Reverse strand** (SubjectStart > SubjectEnd): gene end = SubjectEnd (lower coord)

Positive offsets mean "past the stop codon" for both strands. Window formulas:
- **Forward**: `[gene_end + region_start, gene_end + region_end]`
- **Reverse**: `[gene_end − region_end, gene_end − region_start]` (offsets mirrored so downstream = lower coords)

With `--region_start -50 --region_end 150` both strands get 50 bp overlap into the gene body and 150 bp downstream of the stop codon.

### strand_mode and reading frames
- `same_strand`: only reading frames matching the hit strand (3 forward or 3 reverse)
- `other_strand`: only reading frames on the opposite strand
- `both`: all 6 reading frames (default)

### GTG/TTG start codon handling
`find_orfs.py` records the actual start codon in the ORF DNA. `translate_candidates.py` always translates the first codon as Met (M), regardless of codon identity. This mirrors standard bacterial translation biology (GTG/TTG function as fMet initiators).

### ORF deduplication
`COLLECT_ORFS` merges ORFs from all genome × hit combinations and removes exact DNA duplicates. `TRANSLATE_CANDIDATES` then removes synonymous protein duplicates (different DNA sequences encoding the same protein). Both steps reduce redundancy before CD-HIT clustering.

### CD-HIT parameters
CD-HIT clusters proteins by sequence identity (`--cdhit_identity`, default 0.6) with a minimum alignment length ratio (`--cdhit_length`, default 0.6 applied to both aL and aS). Word size is set to 3 (`-n 3`) to match the 0.6 identity threshold.

---

## Parameters Reference

| Parameter | Default | Description |
|---|---|---|
| `--allele_search_dir` | required | Path to allele_search `blast_summary_filtered/` directory (CSV files) |
| `--genome_dir` | required | Directory with genome FASTA files |
| `--code_dir` | `${projectDir}/code` | Directory with Python scripts |
| `--outdir` | `./results` | Output directory |
| `--region_start` | required | bp offset from gene end for window start (negative = overlaps gene body) |
| `--region_end` | required | bp offset from gene end for window end (negative = overlaps gene body) |
| `--strand_mode` | `both` | `same_strand`, `other_strand`, or `both` |
| `--min_orf_length` | `30` | Minimum ORF length in bp (= 10 aa) |
| `--max_orf_length` | `null` | Maximum ORF length in bp excluding stop codon (optional) |
| `--cdhit_identity` | `0.6` | CD-HIT sequence identity threshold |
| `--cdhit_length` | `0.6` | CD-HIT minimum alignment/length ratio (aL and aS) |

---

## Common Gotchas

- **Strand-aware gene end**: `gene_end` = SubjectEnd for **both** strands (always the stop-codon coordinate). For reverse strand the window offsets are mirrored so positive `region_end` still means "downstream past stop codon". See `extract_downstream.py` for details.
- **Negative `region_start`**: Means the window overlaps the gene body. This is valid (e.g., for overlapping regulatory peptides). `extract_downstream.py` handles this correctly.
- **Contig boundary clipping**: If the region window extends past a contig boundary, the sequence is silently clipped. Short flanking sequences may yield no ORFs — this is expected.
- **ORF deduplication before translation**: `COLLECT_ORFS` deduplicates by exact DNA; `TRANSLATE_CANDIDATES` further deduplicates synonymous protein variants. Both steps are required.
- **GTG/TTG → Met in protein**: `find_orfs.py` outputs the original DNA codon; `translate_candidates.py` always translates the first codon as Met regardless.
- **`--extract` mode requires `--allele_search_dir` and `--genome_dir`**: Needed for region extraction.
- **`-resume` reuses cached tasks**: Omit to force full re-run.
- **Nextflow work directory grows large**: Run `nextflow clean -f` periodically.

---

## Development Workflow

### Adding a new parameter
1. Add to `params {}` in `nextflow.config` with a sensible default
2. Reference as `${params.param_name}` in the relevant module `script:` block
3. Add to the `--help` message in `main.nf`
4. Document in `README.md` and this file's Parameters Reference table

### Adding a new pipeline step
1. Create or edit the relevant `modules/*.nf` file
2. Import the process at the top of `main.nf`
3. Wire into the appropriate workflow block(s) (`full` and any relevant step modes)
4. Update channel routing if downstream steps depend on its output

### Modifying Python scripts
Scripts are called from Nextflow `script:` blocks via:
```
python ${params.code_dir}/script_name.py [args]
```
All scripts use `argparse`. Header parsing/building must go through `header_utils.py`.

### Testing a step change
```bash
# Re-run only the affected step; -resume reuses upstream cached work
nextflow run main.nf --cluster --outdir ./results/test -resume -profile docker
```

---

## Output Interpretation

### Key outputs: `results/cdhit/`

| Path | Description |
|---|---|
| `all_proteins_cdhit.fasta` | Representative protein sequence for each cluster |
| `all_proteins_cdhit.fasta.clstr` | CD-HIT cluster membership file: lists all members per cluster and marks the representative with `*` |

### Intermediate outputs
| Directory | Contents |
|---|---|
| `regions/` | Extracted genomic window FASTA files (one per genome × hit) |
| `candidate_orfs/candidate_orfs.fasta` | Deduplicated candidate ORF DNA sequences |
| `candidate_orfs/candidate_orfs_protein.fasta` | Deduplicated protein sequences (input to CD-HIT) |

### Splitting clusters for downstream use
Use `code/split_clusters.py` to split the merged protein FASTA into one FASTA per cluster, ranked by cluster size:
```bash
python code/split_clusters.py \
  results/cdhit/all_proteins_cdhit.fasta.clstr \
  results/candidate_orfs/candidate_orfs_protein.fasta \
  -o results/clusters/
# Output: cluster_1.fasta (largest), cluster_2.fasta, ...
```

---

## Related Documentation

- `README.md` — User-facing installation and usage guide
- `../allele_search/CLAUDE.md` — Reference for shared conventions
- `TODO.md` — Development task tracking
