"""
Parse a CD-HIT .clstr file and split the original merged protein FASTA into
one FASTA per cluster. Clusters are ranked by size (number of member sequences)
so the largest cluster is named cluster_1.

Output:
  <outdir>/cluster_1.fasta, cluster_2.fasta, ...
"""

import argparse
import re
from pathlib import Path
from Bio import SeqIO


def parse_clstr(clstr_path):
    """Return a list of clusters, each a list of member sequence IDs.
    Order matches the CD-HIT cluster file (cluster 0, 1, 2, ...).
    """
    clusters = []
    current = []
    with open(clstr_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>Cluster'):
                if current:
                    clusters.append(current)
                current = []
            else:
                # e.g.:  0	54aa, >some_id... *
                m = re.search(r'>(.+?)\.\.\.', line)
                if m:
                    current.append(m.group(1))
    if current:
        clusters.append(current)
    return clusters


def main():
    parser = argparse.ArgumentParser(description='Split CD-HIT clusters into per-cluster FASTAs')
    parser.add_argument('clstr_file', help='CD-HIT .clstr file')
    parser.add_argument('fasta_file', help='Merged protein FASTA (representative seqs from cd-hit)')
    parser.add_argument('-o', '--outdir', default='.', help='Output directory')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse clusters
    clusters = parse_clstr(args.clstr_file)
    if not clusters:
        print("WARNING: No clusters found in .clstr file", flush=True)
        return

    # Rank clusters by size descending → cluster_1 is largest
    clusters_ranked = sorted(clusters, key=lambda c: len(c), reverse=True)

    # Build lookup: seq_id → cluster rank (1-based)
    id_to_rank = {}
    for rank, members in enumerate(clusters_ranked, 1):
        for seq_id in members:
            id_to_rank[seq_id] = rank

    # Read merged FASTA and route each sequence to its cluster file.
    # NOTE: cd-hit outputs only the representative sequences in the .fasta file,
    # but the .clstr file contains ALL member IDs. We need the original merged
    # FASTA (all_proteins.fasta) passed as fasta_file, not the cd-hit output.
    handles = {}
    try:
        for rec in SeqIO.parse(args.fasta_file, 'fasta'):
            rank = id_to_rank.get(rec.id)
            if rank is None:
                print(f"WARNING: Sequence '{rec.id}' not found in any cluster — skipping", flush=True)
                continue
            cluster_name = f"cluster_{rank}"
            if cluster_name not in handles:
                out_path = outdir / f"{cluster_name}.fasta"
                handles[cluster_name] = open(out_path, 'w')
            handles[cluster_name].write(f">{rec.description}\n")
            seq = str(rec.seq)
            for i in range(0, len(seq), 60):
                handles[cluster_name].write(seq[i:i+60] + "\n")
    finally:
        for fh in handles.values():
            fh.close()

    n_clusters = len(handles)
    print(f"Split into {n_clusters} cluster(s): " +
          ", ".join(sorted(handles.keys(), key=lambda n: int(n.split('_')[1]))),
          flush=True)


if __name__ == '__main__':
    main()
