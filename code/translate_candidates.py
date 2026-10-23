"""
translate_candidates.py — Translate candidate DNA ORFs to protein and deduplicate.

NEW script for SHP. Called by TRANSLATE_CANDIDATES (exact_match.nf) before the
tblastn exact match search.

Each ORF is translated in frame 0 only — candidate ORFs from find_orfs.py already
start at a valid start codon, so frame 0 is the correct reading frame.  The first
codon is always written as Met (M) regardless of codon identity, matching the
bacterial fMet convention for GTG/TTG initiators (same policy as translate.py).

Identical protein sequences are deduplicated; only the first representative is kept.
This means that multiple DNA ORFs encoding the same protein (synonymous variants)
are collapsed into a single query before tblastn, which then finds ALL genomes
encoding that protein in one search.

Usage:
    python translate_candidates.py \\
        --input  candidate_orfs.fasta \\
        --output candidate_orfs_protein.fasta
"""

import argparse
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from header_utils import parse_header, build_header


def translate_frame0(dna_str):
    """Translate a DNA string in frame 0; return protein without stop codon.

    Strips gap characters and ambiguous Ns, trims to a multiple of 3, then
    translates with BioPython.  The first amino acid is replaced with M to
    handle GTG/TTG start codons (bacterial fMet convention).

    Returns None if the sequence is too short to translate.
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
    # GTG/TTG → Met (first codon is always fMet in bacteria)
    protein = 'M' + protein[1:]
    return protein


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',  required=True,
                        help='candidate_orfs.fasta (DNA, from COLLECT_ORFS)')
    parser.add_argument('--output', required=True,
                        help='candidate_orfs_protein.fasta (protein, for tblastn)')
    args = parser.parse_args()

    seen = set()
    written = 0
    skipped_dup = 0
    skipped_short = 0

    with open(args.output, 'w') as out_fh:
        for record in SeqIO.parse(args.input, 'fasta'):
            protein = translate_frame0(str(record.seq))
            if protein is None:
                skipped_short += 1
                continue
            if protein in seen:
                skipped_dup += 1
                continue
            seen.add(protein)

            seq_id, meta = parse_header(record.description)
            header = build_header(seq_id, meta)
            out_fh.write(header + '\n')
            for i in range(0, len(protein), 60):
                out_fh.write(protein[i:i+60] + '\n')
            written += 1

    print(
        f"INFO [translate_candidates]: {written} unique protein sequences written "
        f"({skipped_dup} duplicate(s) removed, {skipped_short} too-short skipped)",
        flush=True,
    )


if __name__ == '__main__':
    main()
