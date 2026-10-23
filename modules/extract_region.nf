// SHP — extract_region.nf
// Processes: PARSE_ALLELE_INPUT, EXTRACT_DOWNSTREAM_REGION
//
// Step 1: Read allele_search blast_summary_filtered/*.csv → TSV of per-genome gene coordinates
// Step 2: For each genome × hit, extract [gene_end + region_start, gene_end + region_end]
//         from the genome FASTA; strand-aware; clips silently to contig boundaries.

process PARSE_ALLELE_INPUT {
    input:
    path allele_search_dir

    output:
    path "hits.tsv"

    script:
    """
    cat > _parse_allele.py << 'PYEOF'
import csv
from pathlib import Path

def species_from_id(seq_id):
    if 'vs_' in seq_id:
        return seq_id.split('vs_', 1)[1]
    return seq_id

sep = chr(9)   # tab
nl  = chr(10)  # newline
asd = "${allele_search_dir}"

with open('hits.tsv', 'w') as out:
    out.write(sep.join(['source_gene', 'species_name', 'contig', 'subject_start', 'subject_end']) + nl)
    for p in sorted(Path(asd).glob('*.csv')):
        gene = p.stem
        with open(p, newline='') as f:
            for row in csv.DictReader(f):
                sp = species_from_id(row['ID'])
                out.write(sep.join([
                    gene, sp, row['Contig'],
                    row['Subject Start'], row['Subject End']
                ]) + nl)
PYEOF
    python3 _parse_allele.py
    """
}


process EXTRACT_DOWNSTREAM_REGION {
    tag "${species_names[0]}+${species_names.size() - 1}_more"
    publishDir "${params.outdir}/regions", mode: 'copy'

    input:
    tuple path(genome_fastas), val(species_names), val(hits_batch)
    // genome_fastas: list of staged genome files (same index order as species_names, hits_batch)
    // species_names: list of species name strings
    // hits_batch:    list of hits lists, each [[contig, start, end, gene], ...]

    output:
    path "*_region.fasta", optional: true   // multiple files emitted (one per hit across all genomes)

    script:
    // Zip species_names, genome_fastas, hits_batch by index; build all shell commands in one pass
    def cmds = [species_names, genome_fastas, hits_batch]
        .transpose()
        .collect { sp, genome, hits ->
            hits.collect { h ->
                def (contig, start, end_coord, gene) = h
                def out = "${sp}_${gene}_${start}_region.fasta"
                """python ${params.code_dir}/extract_downstream.py \\
        --genome '${genome}' \\
        --contig '${contig}' \\
        --subject_start ${start} \\
        --subject_end ${end_coord} \\
        --region_start ${params.region_start} \\
        --region_end ${params.region_end} \\
        --source_gene '${gene}' \\
        --output '${out}'
    [ -s '${out}' ] || rm -f '${out}'"""
            }.join('\n')
        }.join('\n')
    """
    ${cmds}
    """
}
