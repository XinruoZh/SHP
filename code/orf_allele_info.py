"""
orf_allele_info.py — Translate ORFs, group by protein sequence, and write allele_info output.

Replaces the TRANSLATE_CANDIDATES step. Called by ORF_ALLELE_INFO in analyze_orfs.nf.

Usage:
    python orf_allele_info.py \
        --input  all_orfs_raw.fasta \
        --allele_info_dir allele_info \
        --per_genome_dir  per_genome_orfs

Output (per source_gene):
    allele_info/{source_gene}/self_reference/allele_info_no_strain_name.csv
    allele_info/{source_gene}/self_reference/allele_info_with_strain_name.csv
    allele_info/{source_gene}/self_reference/allele_sequence_output.fasta
    per_genome_orfs/{source_gene}.csv

The allele_info CSV format exactly matches allele_search protein_allele_variant.py output.
Fields that do not apply (variant positions) are written as 'na'.

Allele grouping: unique protein sequences (GTG/TTG start → Met, same as translate_candidates.py).
species_count / species_percentage: unique genome filenames per allele / total unique genomes
for that source_gene (not total ORF count).
strain_names: semicolon-separated sorted genome file basenames (from genome= header field).

Memory design: two-pass streaming.
  Pass 1 — translate each ORF, accumulate only compact protein→genome sets (no DNA retained),
            write orfs_hmmer.fasta inline as HMMER-hit records are encountered.
  Between passes — assign allele names, write allele_info CSVs and FASTA.
  Pass 2 — stream FASTA again, translate, write per-genome CSV rows directly (no accumulation).
This keeps peak RAM proportional to unique proteins × genomes, not to total ORF count.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from header_utils import parse_header


# ---------------------------------------------------------------------------
# Translation (same logic as translate_candidates.py)
# ---------------------------------------------------------------------------

def translate_frame0(dna_str):
    """Translate DNA in frame 0; return protein without stop codon, or None if too short.

    GTG/TTG start codons are translated as Met (bacterial fMet convention).
    """
    dna = dna_str.upper().replace('-', '').replace('N', '')
    dna = dna[:len(dna) - len(dna) % 3]
    if len(dna) < 3:
        return None
    protein = str(Seq(dna).translate())
    if protein.endswith('*'):
        protein = protein[:-1]
    if not protein:
        return None
    protein = 'M' + protein[1:]
    return protein


# ---------------------------------------------------------------------------
# CSV / FASTA writers (matching allele_search format exactly)
# ---------------------------------------------------------------------------

CSV_FIELDNAMES_NO_STRAIN = [
    "allele_name",
    "species_count",
    "species_percentage",
    "variant_position_count",
    "variant_percentage",
    "variant_positions(position:reference\u2192changed)",
    "allele_sequence",
]

CSV_FIELDNAMES_WITH_STRAIN = CSV_FIELDNAMES_NO_STRAIN + ["strain_names"]

HMMER_ALLELE_FIELDS = ["hmmer_hit_count", "hmmer_hit_percentage"]

PER_GENOME_HMMER_FIELDS = ["hmmer_hit", "hmmer_evalue"]


# ---------------------------------------------------------------------------
# HMMER hit table parser
# ---------------------------------------------------------------------------

def parse_hmmer_hits(tblout_path):
    """Parse a merged HMMER --tblout file.

    Returns dict: target_name -> (evalue, score).
    target_name is the full sequence name from the protein FASTA header
    (= record.description in BioPython, i.e. header text without '>').
    Only the best (lowest E-value) hit per target is kept.
    """
    hits = {}
    with open(tblout_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('#') or not line:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            target_name = parts[0]
            try:
                evalue = float(parts[4])
                score  = float(parts[5])
            except ValueError:
                continue
            if target_name not in hits or evalue < hits[target_name][0]:
                hits[target_name] = (evalue, score)
    return hits


# ---------------------------------------------------------------------------
# CSV / FASTA writers (matching allele_search format exactly)
# ---------------------------------------------------------------------------

def _write_csv(allele_list, output_path, with_strain_name, with_hmmer=False):
    """Write allele_info CSV.

    allele_list: list of dicts with keys matching CSV_FIELDNAMES_WITH_STRAIN (+ HMMER fields).
    """
    fieldnames = CSV_FIELDNAMES_WITH_STRAIN if with_strain_name else CSV_FIELDNAMES_NO_STRAIN
    if with_hmmer:
        fieldnames = fieldnames + HMMER_ALLELE_FIELDS
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in allele_list:
            out = {k: row[k] for k in fieldnames}
            writer.writerow(out)


def _write_fasta(allele_list, output_path):
    """Write allele_sequence_output.fasta matching allele_search header format."""
    with open(output_path, "w") as fh:
        for row in allele_list:
            header = (
                f">{row['allele_name']}"
                f"|species_count={row['species_count']}"
                f"|species_percentage={row['species_percentage']}"
            )
            fh.write(header + "\n")
            seq = row["allele_sequence"]
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def _write_fasta_record(fh, header, seq, line_len=60):
    """Write a single FASTA record to an already-open file handle."""
    fh.write(f">{header}\n")
    for i in range(0, len(seq), line_len):
        fh.write(seq[i:i + line_len] + "\n")


# ---------------------------------------------------------------------------
# Per-genome CSV
# ---------------------------------------------------------------------------

PER_GENOME_FIELDNAMES = [
    "genome", "contig",
    "relative_position_start", "relative_position_end",
    "orf_start", "orf_end", "strand",
    "source_gene", "allele_name", "protein_sequence",
]


def _relative_positions(orf_start, orf_end, gene_end, gene_strand):
    """Return (rel_start, rel_end) in the same signed-offset space as --region_start/--region_end.

    Positive = downstream of stop codon. rel_start <= rel_end always.
    Returns ("NA", "NA") when gene_end is None (field absent from header).
    """
    if gene_end is None:
        return "NA", "NA"
    lo, hi = min(orf_start, orf_end), max(orf_start, orf_end)
    if gene_strand == "+":
        return lo - gene_end, hi - gene_end
    else:
        return gene_end - hi, gene_end - lo


# ---------------------------------------------------------------------------
# Allele assignment from compact protein→genome-set structures
# ---------------------------------------------------------------------------

def _assign_alleles(gene_name, protein_to_genomes_local, total_genomes_count,
                    protein_to_hmmer_local, with_hmmer):
    """Group proteins into alleles; return (allele_map, allele_list).

    allele_map: protein → allele_name  (for Pass-2 per-genome CSV lookup)
    allele_list: list of allele row dicts ready for _write_csv / _write_fasta

    protein_to_genomes_local: protein → set of genome strings
    total_genomes_count: int — denominator for species_percentage
    protein_to_hmmer_local: protein → set of genome strings with HMMER hits
                             (empty dict when with_hmmer=False)
    with_hmmer: bool — whether to add hmmer_hit_count / hmmer_hit_percentage columns
    """
    sorted_proteins = sorted(
        protein_to_genomes_local.items(),
        key=lambda x: (-len(x[1]), x[0])
    )

    allele_map = {}
    allele_list = []
    for idx, (protein, genome_set) in enumerate(sorted_proteins, 1):
        allele_name = f"{gene_name}_{idx}"
        allele_map[protein] = allele_name
        unique_genomes = sorted(genome_set)
        species_count = len(unique_genomes)
        species_pct = (species_count / total_genomes_count) * 100 if total_genomes_count else 0.0
        row = {
            "allele_name":   allele_name,
            "species_count": species_count,
            "species_percentage": f"{species_pct:.2f}",
            "variant_position_count": "na",
            "variant_percentage": "na",
            "variant_positions(position:reference\u2192changed)": "",
            "allele_sequence": protein,
            "strain_names": "; ".join(unique_genomes),
        }
        if with_hmmer:
            hit_genomes = protein_to_hmmer_local.get(protein, set())
            hmmer_hit_count = len(hit_genomes)
            row["hmmer_hit_count"] = hmmer_hit_count
            row["hmmer_hit_percentage"] = (
                f"{(hmmer_hit_count / total_genomes_count * 100) if total_genomes_count else 0.0:.2f}"
            )
        allele_list.append(row)

    return allele_map, allele_list


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True,
                        help="Concatenated raw ORF FASTA (all genomes, no dedup)")
    parser.add_argument("--allele_info_dir", default="allele_info",
                        help="Output directory for allele_info format (default: allele_info)")
    parser.add_argument("--per_genome_dir", default="per_genome_orfs",
                        help="Output directory for per-genome CSVs (default: per_genome_orfs)")
    parser.add_argument("--hmmer_hits", default=None,
                        help="Merged HMMER --tblout file (optional)")
    parser.add_argument("--candidate_orfs_dir", default="candidate_orfs",
                        help="Output directory for orfs_hmmer.fasta (default: candidate_orfs)")
    args = parser.parse_args()

    # Load HMMER hits (small dict — stays in memory throughout)
    hmmer_hits = None
    if args.hmmer_hits:
        hmmer_hits = parse_hmmer_hits(args.hmmer_hits)
        print(f"INFO: {len(hmmer_hits)} HMMER hit(s) loaded from {args.hmmer_hits}.", flush=True)

    with_hmmer = hmmer_hits is not None

    # -------------------------------------------------------------------------
    # Pass 1: stream FASTA, translate, build compact gene structures.
    #         Write orfs_hmmer.fasta immediately for HMMER-hit records.
    #         DNA sequence is NOT retained after translation.
    # -------------------------------------------------------------------------
    # gene → protein → set of genome strings (all ORFs)
    protein_to_genomes = defaultdict(lambda: defaultdict(set))
    # gene → protein → set of genome strings where that protein had a HMMER hit
    protein_to_hmmer   = defaultdict(lambda: defaultdict(set))
    # gene → set of all genomes seen (denominator for species_percentage)
    total_genomes_map  = defaultdict(set)

    hmmer_fh = None
    hmmer_fasta_path = Path(args.candidate_orfs_dir) / "orfs_hmmer.fasta"
    if with_hmmer:
        Path(args.candidate_orfs_dir).mkdir(parents=True, exist_ok=True)
        hmmer_fh = open(hmmer_fasta_path, "w")

    total_orfs = 0
    hmmer_hit_count = 0

    for record in SeqIO.parse(args.input, "fasta"):
        seq_id, meta = parse_header(record.description)
        gene     = meta.get("source_gene", "unknown")
        genome   = meta.get("genome", "unknown")
        full_hdr = record.description

        total_orfs += 1
        total_genomes_map[gene].add(genome)

        protein = translate_frame0(str(record.seq))   # DNA discarded after this line
        if protein is not None:
            protein_to_genomes[gene][protein].add(genome)
            if with_hmmer and full_hdr in hmmer_hits:
                protein_to_hmmer[gene][protein].add(genome)

        # Write HMMER-hit DNA to orfs_hmmer.fasta immediately (no buffering)
        if with_hmmer and full_hdr in hmmer_hits:
            _write_fasta_record(hmmer_fh, full_hdr, str(record.seq))
            hmmer_hit_count += 1

    if hmmer_fh:
        hmmer_fh.close()
        print(f"INFO: {hmmer_hit_count} HMMER-hit ORF(s) written to {hmmer_fasta_path}",
              flush=True)

    print(f"INFO: {total_orfs} ORF(s) read across {len(total_genomes_map)} source gene(s).",
          flush=True)

    # -------------------------------------------------------------------------
    # Between passes: assign allele names, write allele_info CSVs and FASTA.
    # -------------------------------------------------------------------------
    allele_maps = {}  # gene → protein → allele_name  (kept for Pass 2 lookup)

    for gene_name in sorted(total_genomes_map.keys()):
        p2g = protein_to_genomes[gene_name]
        if not p2g:
            print(f"WARNING [{gene_name}]: No translatable ORFs — skipping output.", flush=True)
            continue

        total_n = len(total_genomes_map[gene_name])
        p2h = protein_to_hmmer[gene_name]  # empty defaultdict(set) when no HMMER hits

        allele_map, allele_list = _assign_alleles(
            gene_name, p2g, total_n, p2h, with_hmmer
        )
        allele_maps[gene_name] = allele_map

        # Write self_reference outputs
        out_dir = Path(args.allele_info_dir) / gene_name / "self_reference"
        out_dir.mkdir(parents=True, exist_ok=True)
        _write_csv(allele_list, out_dir / "allele_info_no_strain_name.csv",
                   with_strain_name=False, with_hmmer=with_hmmer)
        _write_csv(allele_list, out_dir / "allele_info_with_strain_name.csv",
                   with_strain_name=True, with_hmmer=with_hmmer)
        _write_fasta(allele_list, out_dir / "allele_sequence_output.fasta")
        print(f"INFO [{gene_name}]: {len(allele_list)} allele(s) written to {out_dir}",
              flush=True)

        # Write hmmer_only_self_reference (proteins where ≥1 genome had a HMMER hit)
        if with_hmmer and p2h:
            _, hmmer_allele_list = _assign_alleles(
                gene_name, p2h, total_n, p2h, with_hmmer=True
            )
            if hmmer_allele_list:
                hmmer_dir = Path(args.allele_info_dir) / gene_name / "hmmer_only_self_reference"
                hmmer_dir.mkdir(parents=True, exist_ok=True)
                _write_csv(hmmer_allele_list, hmmer_dir / "allele_info_no_strain_name.csv",
                           with_strain_name=False, with_hmmer=True)
                _write_csv(hmmer_allele_list, hmmer_dir / "allele_info_with_strain_name.csv",
                           with_strain_name=True, with_hmmer=True)
                _write_fasta(hmmer_allele_list, hmmer_dir / "allele_sequence_output.fasta")
                print(f"INFO [{gene_name}]: {len(hmmer_allele_list)} HMMER-hit allele(s) "
                      f"written to {hmmer_dir}", flush=True)

    if not allele_maps:
        print("WARNING: No source genes had translatable ORFs.", flush=True)
        print("INFO: orf_allele_info complete.", flush=True)
        return

    # -------------------------------------------------------------------------
    # Pass 2: stream FASTA again, write per-genome CSV rows directly.
    #         One csv.DictWriter per gene is kept open; rows written immediately.
    # -------------------------------------------------------------------------
    Path(args.per_genome_dir).mkdir(parents=True, exist_ok=True)
    per_genome_fieldnames = (
        PER_GENOME_FIELDNAMES + PER_GENOME_HMMER_FIELDS if with_hmmer
        else PER_GENOME_FIELDNAMES
    )

    per_genome_fhs      = {}
    per_genome_writers  = {}
    for gene_name in sorted(allele_maps.keys()):
        fh = open(Path(args.per_genome_dir) / f"{gene_name}.csv", "w", newline="")
        writer = csv.DictWriter(fh, fieldnames=per_genome_fieldnames)
        writer.writeheader()
        per_genome_fhs[gene_name]     = fh
        per_genome_writers[gene_name] = writer

    for record in SeqIO.parse(args.input, "fasta"):
        seq_id, meta = parse_header(record.description)
        gene     = meta.get("source_gene", "unknown")
        genome   = meta.get("genome", "unknown")
        full_hdr = record.description

        if gene not in per_genome_writers:
            continue

        protein = translate_frame0(str(record.seq))
        if protein is None:
            continue

        allele_name = allele_maps[gene].get(protein, "")

        try:
            orf_start = int(meta["SubjectStart"])
            orf_end   = int(meta["SubjectEnd"])
        except (KeyError, ValueError):
            orf_start = orf_end = 0

        orf_strand  = "+" if orf_start <= orf_end else "-"
        gene_end_str = meta.get("gene_end")
        gene_end   = int(gene_end_str) if gene_end_str is not None else None
        gene_strand = meta.get("gene_strand", "+")
        strand      = "+" if orf_strand == gene_strand else "-"
        contig     = meta.get("contig", "unknown")

        rel_start, rel_end = _relative_positions(orf_start, orf_end, gene_end, gene_strand)

        row = {
            "genome":                  genome,
            "contig":                  contig,
            "relative_position_start": rel_start,
            "relative_position_end":   rel_end,
            "orf_start":               orf_start,
            "orf_end":                 orf_end,
            "strand":                  strand,
            "source_gene":             gene,
            "allele_name":             allele_name,
            "protein_sequence":        protein,
        }
        if with_hmmer:
            is_hit = full_hdr in hmmer_hits
            row["hmmer_hit"]    = is_hit
            row["hmmer_evalue"] = hmmer_hits[full_hdr][0] if is_hit else "NA"

        per_genome_writers[gene].writerow(row)

    for fh in per_genome_fhs.values():
        fh.close()

    for gene_name in sorted(allele_maps.keys()):
        csv_path = Path(args.per_genome_dir) / f"{gene_name}.csv"
        print(f"INFO [{gene_name}]: per-genome CSV written to {csv_path}", flush=True)

    print("INFO: orf_allele_info complete.", flush=True)


if __name__ == "__main__":
    main()
