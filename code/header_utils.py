"""
header_utils.py

Canonical parser and builder for the internal pipeline FASTA header format:

    >ID|contig=NAME|E-value=VALUE|AlignLen=LEN|SeqIdentity=%|SubjectStart=POS|SubjectEnd=POS[|source_gene=GENE]

The optional ``source_gene`` field is new in SHP (absent in allele_search).
It tracks which allele_search target gene produced each ORF candidate.

This format is written by csv_to_fasta.py and read by start_stop_codon_dna_seq.py,
up_down_stream.py, extract_dna_from_blast_csv.py, and protein_allele_variant.py.
All header parsing and building should go through these functions so that a format
change only requires editing this file.
"""

# Canonical field order — all writers follow this so readers can rely on it.
# source_gene is optional (SHP-only); build_header emits it only when present.
HEADER_FIELDS = ('contig', 'E-value', 'AlignLen', 'SeqIdentity', 'SubjectStart', 'SubjectEnd', 'source_gene')


def parse_header(header: str) -> tuple:
    """Parse a pipeline FASTA header into (seq_id, metadata_dict).

    Handles an optional leading '>'.  Each field is whitespace-stripped.
    Parts that do not contain '=' are skipped with a printed warning so
    that callers never receive partial data silently.

    All metadata values are returned as strings; callers cast as needed
    (e.g. int(metadata['SubjectStart'])).

    The ``source_gene`` field is optional; it will be present in the returned
    dict only when the header contains it.

    Returns
    -------
    seq_id   : str  — the ID before the first '|'
    metadata : dict — key/value pairs from the remaining '|'-separated fields
    """
    header = header.lstrip('>').strip()
    parts = header.split('|')
    seq_id = parts[0].strip()
    metadata = {}
    for part in parts[1:]:
        part = part.strip()
        if not part:
            continue
        if '=' in part:
            key, value = part.split('=', 1)
            metadata[key.strip()] = value.strip()
        else:
            print(f"WARNING: Ignoring malformed header field '{part}' in: {header}", flush=True)
    return seq_id, metadata


def build_header(seq_id: str, metadata: dict) -> str:
    """Rebuild a pipeline FASTA header string from seq_id and a metadata dict.

    Fields are emitted in canonical order (HEADER_FIELDS).  Any keys not in
    HEADER_FIELDS are appended afterwards so no data is silently dropped.

    The ``source_gene`` field is only emitted when present in metadata.

    The returned string always begins with '>'.
    """
    parts = [f">{seq_id}"]
    for key in HEADER_FIELDS:
        if key in metadata:
            parts.append(f"{key}={metadata[key]}")
    for key, value in metadata.items():
        if key not in HEADER_FIELDS:
            parts.append(f"{key}={value}")
    return '|'.join(parts)


def species_from_id(seq_id: str) -> str:
    """Return the genome/species name embedded in a BLAST sequence ID.

    BLAST output files are named '{query}_vs_{genome}', so every sequence ID
    produced by the pipeline carries that same prefix.  This function strips
    everything up to and including the first 'vs_' token.

    A leading '>' is removed if present.  If no 'vs_' token is found the
    original ID is returned unchanged (e.g. when the input is already a
    bare genome name).

    Examples
    --------
    >>> species_from_id('insertion_region_1_vs_strain_ABC')
    'strain_ABC'
    >>> species_from_id('>query_vs_genome_vs_suffix')
    'genome_vs_suffix'          # preserves everything after the first vs_
    >>> species_from_id('plain_genome_name')
    'plain_genome_name'
    """
    seq_id = seq_id.lstrip('>').strip()
    if 'vs_' in seq_id:
        return seq_id.split('vs_', 1)[1]
    return seq_id
