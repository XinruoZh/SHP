"""
extract_downstream.py — Extract a genomic window relative to an allele_search gene hit.

NEW script (not in allele_search). Called by EXTRACT_DOWNSTREAM_REGION in extract_region.nf.

Usage:
    python extract_downstream.py \
        --genome <genome.fasta> \
        --contig <contig_name> \
        --subject_start <int> \
        --subject_end <int> \
        --region_start <int> \
        --region_end <int> \
        --source_gene <query_name> \
        --output <out.fasta>

Logic:
    Strand detection: SubjectStart > SubjectEnd → reverse strand (same convention as allele_search).

    Gene end (strand-aware):
        Forward strand (SubjectStart <= SubjectEnd): gene_end = SubjectEnd
        Reverse strand (SubjectStart >  SubjectEnd): gene_end = SubjectStart  ← the larger coord

    Window: [gene_end + region_start, gene_end + region_end]
        Both offsets may be negative (window overlaps gene body).
        Window is clipped to [0, contig_length] if out of bounds (no error, log to stderr).

    Always extracts the raw forward-strand genomic sequence. The gene_strand field in
    the output header encodes the original hit strand; find_orfs.py uses it together
    with --strand_mode to select which reading frames to scan.
"""

import argparse
import os
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from header_utils import build_header


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True,
                        help="Genome FASTA file")
    parser.add_argument("--contig", required=True,
                        help="Contig/chromosome name in that genome")
    parser.add_argument("--subject_start", type=int, required=True,
                        help="BLAST SubjectStart (1-based)")
    parser.add_argument("--subject_end", type=int, required=True,
                        help="BLAST SubjectEnd (1-based)")
    parser.add_argument("--region_start", type=int, required=True,
                        help="bp offset from gene end for window start (can be negative)")
    parser.add_argument("--region_end", type=int, required=True,
                        help="bp offset from gene end for window end (can be negative)")
    parser.add_argument("--source_gene", required=True,
                        help="Query gene name from allele_search (e.g. 'murM')")
    parser.add_argument("--output", required=True,
                        help="Output FASTA file path")
    args = parser.parse_args()

    # Load genome into dict keyed by contig ID
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    if args.contig not in genome:
        print(
            f"WARNING: contig '{args.contig}' not found in {args.genome}; skipping.",
            file=sys.stderr,
        )
        # Write empty output so Nextflow does not fail on a missing file
        open(args.output, "w").close()
        return

    contig_seq = genome[args.contig].seq
    contig_len = len(contig_seq)

    # Strand detection: SubjectStart > SubjectEnd → reverse strand
    is_reverse = args.subject_start > args.subject_end

    # Gene end = stop codon position (3' end of transcript), strand-aware:
    #   Forward: SubjectEnd (higher coord) is the stop codon
    #   Reverse: SubjectEnd (lower coord) is the stop codon
    # Window offsets are relative to this anchor.  Positive region_end = downstream
    # past the stop codon.  Both offsets may be negative (overlap with gene body).
    #
    # For reverse-strand genes, "downstream" is at LOWER genomic coordinates,
    # so the offsets are mirrored: win = [gene_end - region_end, gene_end - region_start].
    gene_end = args.subject_end  # stop codon for both strands

    if is_reverse:
        # Downstream is toward lower coords; flip the signed offsets so that
        # a positive region_end still means "past the stop codon".
        win_start = gene_end - args.region_end
        win_end   = gene_end - args.region_start
    else:
        win_start = gene_end + args.region_start
        win_end   = gene_end + args.region_end

    # Convert to 0-based Python slice; clip to contig boundaries
    slice_start = max(0, win_start - 1)
    slice_end   = min(contig_len, win_end)

    if slice_start >= slice_end:
        print(
            f"WARNING: window [{win_start}, {win_end}] for contig '{args.contig}' "
            f"is empty after clipping to [1, {contig_len}]; skipping.",
            file=sys.stderr,
        )
        open(args.output, "w").close()
        return

    if slice_start > win_start - 1 or slice_end < win_end:
        print(
            f"WARNING: window [{win_start}, {win_end}] clipped to "
            f"[{slice_start + 1}, {slice_end}] for contig '{args.contig}'.",
            file=sys.stderr,
        )

    region_seq = str(contig_seq[slice_start:slice_end])

    # seq_id encodes 1-based genome coords for traceability
    seq_id = f"{args.source_gene}_region_{args.contig}_{slice_start + 1}_{slice_end}"
    genome_name = Path(args.genome).stem

    header = build_header(seq_id, {
        "contig":       args.contig,
        "SubjectStart": slice_start + 1,   # 1-based genome coord of region start
        "SubjectEnd":   slice_end,          # 1-based genome coord of region end
        "gene_strand":  "-" if is_reverse else "+",
        "gene_end":     gene_end,           # 1-based stop-codon coord; anchor for relative positions
        "source_gene":  args.source_gene,
        "genome":       genome_name,
    })

    with open(args.output, "w") as fh:
        fh.write(header + "\n")
        for i in range(0, len(region_seq), 60):
            fh.write(region_seq[i:i+60] + "\n")


if __name__ == "__main__":
    main()
