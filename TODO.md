# SHP Development TODO

## Phase 1: Core scripts

- [x] `code/extract_downstream.py` — Extract genomic window from genome FASTA
      - Read allele_search blast_summary_filtered CSV
      - Compute gene_end (strand-aware: forward=SubjectEnd, reverse=SubjectEnd of lower coord)
      - Extract [gene_end + region_start, gene_end + region_end]
      - Clip to contig boundaries (no error, log to stderr)
      - Output FASTA with header_utils.py header format + source_gene field

- [x] `code/find_orfs.py` — Find ORFs in extracted region
      - Scan reading frames per strand_mode (same_strand/other_strand/both)
      - Start codons: ATG, GTG, TTG; Stop codons: TAA, TAG, TGA
      - Filter by min_orf_length (default 30 bp)
      - Output ORF FASTA (DNA, original start codon)

- [x] `code/translate_candidates.py` — Translate DNA ORFs to protein
      - Frame 0 only (ORFs already start at valid start codon)
      - First codon always written as Met (GTG/TTG → fMet)
      - Deduplicate synonymous protein variants

- [x] `code/split_clusters.py` — Standalone: split CD-HIT .clstr into per-cluster FASTAs
      - Rank clusters by size (largest = cluster_1)
      - Output: cluster_1.fasta, cluster_2.fasta, ...

- [x] `code/header_utils.py` (+ SHP: added `source_gene` to HEADER_FIELDS)

## Phase 2: Nextflow modules

- [x] `modules/extract_region.nf` — PARSE_ALLELE_INPUT, EXTRACT_DOWNSTREAM_REGION
- [x] `modules/orf_find.nf` — FIND_ORFS, COLLECT_ORFS (+ deduplication logic)
- [x] `modules/translate.nf` — TRANSLATE_CANDIDATES
- [x] `modules/cluster.nf` — CD_HIT_CLUSTER

## Phase 3: Workflow wiring

- [x] `main.nf` — Full workflow + step modes (extract/cluster) + validation
- [x] `nextflow.config` — Parameters and resource defaults

## Phase 4: Docker

- [x] `docker/Dockerfile` — Pinned versions; build + smoke-test still needed
- [x] `environment.yml` — Conda fallback environment

## Phase 5: Testing

- [x] Add test fixtures to `test/` — `test/genomes/` (20 symlinks) + `test/README.md`
- [ ] End-to-end test with allele_search example output (Rgg144 or murM)
      - Smoke test: `test/genomes/` (20 strains), `--region_start -200 --region_end 200 --strand_mode other_strand`, outdir `test/results/smoke/`
      - Full test: MitisCollectionPubMLST_Fasta (~500+ genomes), outdir `test/results/full/`
- [ ] Verify CD-HIT output: check that cluster sizes and representatives are sensible
