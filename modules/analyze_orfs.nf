// SHP — analyze_orfs.nf
// Process: ORF_ALLELE_INFO
//
// Step 5: Translate raw ORFs, group by protein sequence per source_gene, and write:
//   - allele_info/{source_gene}/self_reference/  (exact allele_search allele_info format)
//   - per_genome_orfs/{source_gene}.csv          (one row per ORF occurrence with genome info)

process ORF_ALLELE_INFO {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path raw_orfs
    path hmmer_hits  // all_hits.tblout when --hmmer_query is set; [] (NO_FILE) otherwise

    output:
    path "allele_info/**",                          emit: allele_info
    path "per_genome_orfs/**",                      emit: per_genome
    path "candidate_orfs/orfs_hmmer.fasta",         optional: true

    script:
    def hmmer_arg = hmmer_hits.name != 'NO_FILE' ? "--hmmer_hits ${hmmer_hits}" : ""
    """
    python ${params.code_dir}/orf_allele_info.py \\
        --input "${raw_orfs}" \\
        --allele_info_dir allele_info \\
        --per_genome_dir  per_genome_orfs \\
        --candidate_orfs_dir candidate_orfs \\
        ${hmmer_arg}
    """
}
