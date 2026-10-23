// SHP — orf_find.nf
// Processes: FIND_ORFS, COLLECT_ORFS
//
// Step 3: Find all ORFs >= min_orf_length in each extracted region.
//         Reads frames per strand_mode; start codons set by --start_codons (default: ATG).
// Step 4: Concatenate ORFs from all genomes; preserve full headers (genome info retained).

process FIND_ORFS {

    input:
    path region_fastas   // list of up to batch_size region FASTAs, staged in work dir

    output:
    path "batch_orfs.fasta", optional: true

    script:
    """
    touch batch_orfs.fasta
    for f in ${region_fastas}; do
        stem="\${f%_region.fasta}"
        python ${params.code_dir}/find_orfs.py \\
            --input "\$f" \\
            --output "\${stem}_orfs.fasta" \\
            --min_orf_length ${params.min_orf_length} \\
            ${params.max_orf_length != null ? "--max_orf_length ${params.max_orf_length}" : ""} \\
            --strand_mode ${params.strand_mode} \\
            --start_codons ${params.start_codons}
        [ -s "\${stem}_orfs.fasta" ] && cat "\${stem}_orfs.fasta" >> batch_orfs.fasta || true
    done
    [ -s batch_orfs.fasta ] || rm -f batch_orfs.fasta
    """
}


process COLLECT_ORFS {
    publishDir "${params.outdir}/candidate_orfs", mode: 'copy'

    input:
    path orf_fastas, stageAs: "orfs_?.fasta"

    output:
    path "all_orfs_raw.fasta"

    script:
    """
    # Concatenate all ORF FASTAs; preserve full headers (genome info retained for allele_info step)
    for f in orfs_*.fasta; do cat "\$f"; done > all_orfs_raw.fasta
    """
}
