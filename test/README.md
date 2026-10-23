# SHP Test Suite — Rgg144 in S. mitis

## Test Data

| Item | Path |
|---|---|
| allele_search input (237 Rgg144 hits) | `/home/xinruoz/Document/Hiller/allele_search/results/mitis_Rgg144/blast_summary_filtered/` |
| Smoke-test genomes (20 strains) | `test/genomes/` (symlinks into MitisCollectionPubMLST_Fasta/) |
| Full genome collection | `/home/xinruoz/Document/Hiller/mitis/data/MitisCollectionPubMLST_Fasta/` |

`test/genomes/` was created with:
```bash
python3 - <<'EOF'
import csv, os, pathlib
csv_path = "/home/xinruoz/Document/Hiller/allele_search/results/mitis_Rgg144/blast_summary_filtered/Rgg144.csv"
genome_src = "/home/xinruoz/Document/Hiller/mitis/data/MitisCollectionPubMLST_Fasta"
out_dir = pathlib.Path("test/genomes"); out_dir.mkdir(exist_ok=True)
seen, count = set(), 0
with open(csv_path) as f:
    for row in csv.DictReader(f):
        if count >= 20: break
        species = row["ID"].split("vs_", 1)[1]
        if species in seen: continue
        seen.add(species)
        src = pathlib.Path(genome_src) / f"{species}.fasta"
        dst = out_dir / f"{species}.fasta"
        if src.exists() and not dst.exists(): os.symlink(src, dst); count += 1
EOF
```

## Requirements

- Nextflow installed on host
- Docker image built: `docker build -t shp:latest docker/`

## Test Parameters

| Parameter | Value | Rationale |
|---|---|---|
| `--region_start` | `-200` | 200 bp overlap into gene body (3' end of Rgg144) |
| `--region_end` | `200` | 200 bp beyond gene end |
| `--strand_mode` | `other_strand` | sORFs on antisense strand relative to Rgg144 |

## Smoke Test (20 genomes, ~5-10 min)

Fast iterative test. All hits will have `has_target_gene = True` (subset only contains Rgg144 strains).

```bash
cd /home/xinruoz/Document/Hiller/SHP

nextflow run main.nf \
  --allele_search_dir /home/xinruoz/Document/Hiller/allele_search/results/mitis_Rgg144/blast_summary_filtered \
  --genome_dir test/genomes \
  --region_start -200 \
  --region_end 200 \
  --strand_mode other_strand \
  --outdir test/results/smoke \
  -profile docker \
  -resume
```

## Full Integration Test (~500+ genomes)

Searches across the entire mitis collection. Strains without Rgg144 will appear as
`has_target_gene = False` in allele_info CSVs — this verifies the annotation step.

```bash
cd /home/xinruoz/Document/Hiller/SHP

nextflow run main.nf \
  --allele_search_dir /home/xinruoz/Document/Hiller/allele_search/results/mitis_Rgg144/blast_summary_filtered \
  --genome_dir /home/xinruoz/Document/Hiller/mitis/data/MitisCollectionPubMLST_Fasta \
  --region_start -200 \
  --region_end 200 \
  --strand_mode other_strand \
  --outdir test/results/full \
  -profile docker \
  -resume
```

## Verification Checklist

Run after each test. Replace `test/results/smoke` with `test/results/full` for the full test.

### 1. Region extraction
```bash
ls test/results/smoke/regions/ | wc -l        # should equal number of CSV hits processed
head -2 test/results/smoke/regions/*.fasta | head -10
# Check: headers contain source_gene=Rgg144, sequences ~400 bp
```

### 2. ORF candidates
```bash
grep -c "^>" test/results/smoke/candidate_orfs/candidate_orfs.fasta
# Check: non-zero count; headers contain source_gene=Rgg144
grep "^>" test/results/smoke/candidate_orfs/candidate_orfs.fasta | head -5
```

### 3. Exact match results
```bash
ls test/results/smoke/exact_match_filtered/
# Check: at least one non-empty CSV (pipeline skips ORFs with zero hits)
head -2 test/results/smoke/exact_match_filtered/*.csv 2>/dev/null | head -20
```

### 4. Translated proteins
```bash
ls test/results/smoke/protein_seq/inframe/
grep -c "^>" test/results/smoke/protein_seq/inframe/*.fasta 2>/dev/null
```

### 5. Allele info (key output)
```bash
ls test/results/smoke/allele_info/
# Check: subdirectory per ORF name
head -3 test/results/smoke/allele_info/*/consensus_reference/allele_info_no_strain_name.csv
# Check: has_target_gene column present; allele_1 row exists
```

### 6. Annotation check (full test only)
```bash
python3 -c "
import pandas as pd, glob
for f in sorted(glob.glob('test/results/full/allele_info/*/consensus_reference/allele_info_no_strain_name.csv')):
    df = pd.read_csv(f)
    print(f.split('allele_info/')[1].split('/')[0], df['has_target_gene'].value_counts().to_dict())
"
# Expected: mix of True/False per ORF (sORF found in strains without Rgg144)
```

### 7. Alignment figures and logos
```bash
ls test/results/smoke/alignment_figures/*/consensus_reference/*.pdf 2>/dev/null
ls test/results/smoke/logo/ 2>/dev/null
```

## Step-mode Re-runs

After a full run, re-run individual stages without reprocessing upstream steps:

```bash
# Re-run only allele analysis + figures (e.g. after changing species_threshold)
nextflow run main.nf --visual \
  --allele_search_dir /home/xinruoz/Document/Hiller/allele_search/results/mitis_Rgg144/blast_summary_filtered \
  --outdir test/results/smoke -profile docker -resume

# Re-run only annotation
nextflow run main.nf --annotate \
  --allele_search_dir /home/xinruoz/Document/Hiller/allele_search/results/mitis_Rgg144/blast_summary_filtered \
  --outdir test/results/smoke -profile docker -resume
```
