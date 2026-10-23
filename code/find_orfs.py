"""
find_orfs.py — Find all ORFs in an extracted genomic region.

NEW script (not in allele_search). Called by FIND_ORFS in orf_find.nf.

Usage:
    python find_orfs.py \
        --input <region.fasta> \
        --min_orf_length <int>  (default: 30 bp = 10 aa) \
        --strand_mode <same_strand|other_strand|both> \
        --output <orfs.fasta>

Logic:
    Input: FASTA file with one or more region sequences (output of extract_downstream.py).
           Each record header encodes contig, genomic coords, gene_strand, source_gene
           (via header_utils.py format).

    Reading frames scanned are determined by strand_mode + gene_strand from the header:
        same_strand:  frames on the same strand as the allele_search hit
        other_strand: frames on the opposite strand
        both:         all 6 frames (default)

    ORF finding uses a stop-to-stop codon scan across all 3 reading frames in the
    requested strands (Bio.Seq.reverse_complement for the reverse strand).
    Start codons accepted: ATG, GTG, TTG (bacterial fMet initiators).
    Stop codons: TAA, TAG, TGA.

    One ORF is reported per stop codon. If an ATG exists between the previous
    stop and the current stop, the first ATG is used as the start (canonical
    bacterial start). Otherwise the first alternative start codon (GTG/TTG) is used.

    Output sequences are ORF DNA in 5'→3' coding-strand orientation (reverse-strand ORFs
    are output as their reverse complement). The original start codon is kept as-is;
    translate.py always translates the first codon as Met (M) downstream.

    Genome coordinates in the output header use the allele_search convention:
        SubjectStart > SubjectEnd  →  reverse strand
        SubjectStart ≤ SubjectEnd  →  forward strand

    This script does NOT deduplicate. COLLECT_ORFS (Nextflow) deduplicates across all
    genomes before BLAST.
"""

import argparse
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from header_utils import parse_header, build_header

STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
ALL_START_CODONS = {"ATG", "GTG", "TTG"}


def _find_orfs_in_seq(seq, min_len, max_len, start_codons):
    """Scan all 3 reading frames of seq for ORFs.

    One ORF per stop codon. Start codon selection between two consecutive stops:
    - If any ATG is present, use the first ATG (canonical bacterial start).
    - Otherwise, use the first alternative start codon (GTG or TTG).

    Returns list of (s_pos, e_pos, start_codon, orf_dna):
        s_pos       — 0-based start of the start codon
        e_pos       — 0-based exclusive end (last nt of stop codon + 1)
        start_codon — 'ATG', 'GTG', or 'TTG'
        orf_dna     — DNA from start codon through end of stop codon
    """
    orfs = []
    n = len(seq)
    for frame in range(3):
        first_start = None      # first start codon (any) since last stop
        first_atg   = None      # first ATG since last stop
        pos = frame
        while pos + 2 < n:
            codon = seq[pos:pos + 3]
            if codon in start_codons:
                if first_start is None:
                    first_start = (pos, codon)
                if codon == "ATG" and first_atg is None:
                    first_atg = (pos, codon)
            elif codon in STOP_CODONS:
                chosen = first_atg if first_atg is not None else first_start
                if chosen is not None:
                    s_pos, s_codon = chosen
                    coding_len = pos - s_pos    # bp excluding stop codon
                    if coding_len >= min_len and (max_len is None or coding_len <= max_len):
                        e_pos = pos + 3
                        orfs.append((s_pos, e_pos, s_codon, seq[s_pos:e_pos]))
                first_start = None
                first_atg   = None
            pos += 3
    return orfs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",  required=True,
                        help="FASTA from extract_downstream.py")
    parser.add_argument("--min_orf_length", type=int, default=30,
                        help="Minimum ORF length in bp, excluding stop codon (default: 30)")
    parser.add_argument("--max_orf_length", type=int, default=None,
                        help="Maximum ORF length in bp, excluding stop codon (optional)")
    parser.add_argument("--strand_mode", default="both",
                        choices=["same_strand", "other_strand", "both"],
                        help="Which strands to scan (default: both)")
    parser.add_argument("--start_codons", default="ATG",
                        help="Comma-separated start codons to accept (default: ATG)")
    parser.add_argument("--output", required=True,
                        help="Output ORF FASTA file path")
    args = parser.parse_args()

    start_codons = frozenset(c.strip().upper() for c in args.start_codons.split(","))
    invalid = start_codons - ALL_START_CODONS
    if invalid:
        sys.exit(f"ERROR: unrecognised start codon(s): {', '.join(sorted(invalid))}. "
                 f"Allowed: {', '.join(sorted(ALL_START_CODONS))}")

    with open(args.output, "w") as out_fh:
        for record in SeqIO.parse(args.input, "fasta"):
            seq_id, meta = parse_header(record.description)
            contig      = meta.get("contig", "unknown")
            win_start   = int(meta["SubjectStart"])   # 1-based genome coord of region start
            gene_strand = meta.get("gene_strand", "+")
            gene_end    = meta.get("gene_end", "0")   # stop-codon coord; propagated to ORF headers
            source_gene = meta.get("source_gene", "unknown")
            genome      = meta.get("genome", "unknown")

            region_seq = str(record.seq).upper()
            n = len(region_seq)

            # Determine which strands to scan
            if args.strand_mode == "both":
                scan_forward = scan_reverse = True
            elif args.strand_mode == "same_strand":
                scan_forward = (gene_strand == "+")
                scan_reverse = (gene_strand == "-")
            else:   # other_strand
                scan_forward = (gene_strand == "-")
                scan_reverse = (gene_strand == "+")

            # --- Forward strand ---
            if scan_forward:
                for s_pos, e_pos, s_codon, orf_dna in _find_orfs_in_seq(
                        region_seq, args.min_orf_length, args.max_orf_length, start_codons):
                    g_start = win_start + s_pos       # 1-based, 5' end
                    g_end   = win_start + e_pos - 1   # 1-based, 3' end (last nt of stop)
                    # SubjectStart ≤ SubjectEnd  →  forward strand
                    orf_id = f"{source_gene}_{contig}_{g_start}_{g_end}"
                    out_fh.write(build_header(orf_id, {
                        "contig":       contig,
                        "SubjectStart": g_start,
                        "SubjectEnd":   g_end,
                        "source_gene":  source_gene,
                        "genome":       genome,
                        "gene_strand":  gene_strand,
                        "gene_end":     gene_end,
                    }) + "\n")
                    for i in range(0, len(orf_dna), 60):
                        out_fh.write(orf_dna[i:i + 60] + "\n")

            # --- Reverse strand ---
            if scan_reverse:
                rc_seq = str(Seq(region_seq).reverse_complement())
                for s_pos, e_pos, s_codon, orf_dna in _find_orfs_in_seq(
                        rc_seq, args.min_orf_length, args.max_orf_length, start_codons):
                    # Map RC 0-based positions back to 1-based genome coordinates.
                    # RC position p  →  original 0-based position (n - 1 - p).
                    # 5' end of ORF (start codon) = LARGER genome coord (SubjectStart).
                    g_start = win_start + (n - 1 - s_pos)   # 5' end, larger
                    g_end   = win_start + (n - e_pos)        # 3' end, smaller
                    # SubjectStart > SubjectEnd  →  reverse strand
                    orf_id = f"{source_gene}_{contig}_{g_start}_{g_end}"
                    out_fh.write(build_header(orf_id, {
                        "contig":       contig,
                        "SubjectStart": g_start,
                        "SubjectEnd":   g_end,
                        "source_gene":  source_gene,
                        "genome":       genome,
                        "gene_strand":  gene_strand,
                        "gene_end":     gene_end,
                    }) + "\n")
                    for i in range(0, len(orf_dna), 60):
                        out_fh.write(orf_dna[i:i + 60] + "\n")


if __name__ == "__main__":
    main()
