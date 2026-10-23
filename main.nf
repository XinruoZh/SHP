#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// SHP — Short Peptide Hunting Pipeline
// main.nf — Workflow entrypoint
//
// Finds short ORFs (sORFs) in genomic regions flanking allele_search hits,
// translates them to protein, groups by protein sequence, and writes allele_info output.

// ---------------------------------------------------------------------------
// Module imports
// ---------------------------------------------------------------------------
include { PARSE_ALLELE_INPUT; EXTRACT_DOWNSTREAM_REGION } from './modules/extract_region'
include { FIND_ORFS; COLLECT_ORFS }                       from './modules/orf_find'
include { HMMER_SEARCH }                                  from './modules/hmmer_search'
include { ORF_ALLELE_INFO }                               from './modules/analyze_orfs'

// ---------------------------------------------------------------------------
// Help message
// ---------------------------------------------------------------------------
def helpMessage() {
    log.info """
    ╔═══════════════════════════════════════════════════════════╗
    ║   SHP — Short Peptide Hunting Pipeline                    ║
    ╚═══════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf [options] -profile docker

    Required parameters:
        --allele_search_dir   Path to allele_search blast_summary_filtered/ directory
        --genome_dir          Directory with genome FASTA files
        --region_start        bp offset from gene end for window start (can be negative)
        --region_end          bp offset from gene end for window end   (can be negative)

    Optional parameters:
        --strand_mode         same_strand | other_strand | both  [default: both]
        --min_orf_length      Minimum ORF length in bp           [default: 30]
        --max_orf_length      Maximum ORF length in bp, excl. stop codon [optional]
        --batch_size          Region FASTAs per FIND_ORFS task   [default: 100]
        --outdir              Output directory                    [default: ./results]

    HMMER search (step 3.5, optional — skipped if --hmmer_query is not set):
        --hmmer_query         Path to query (.faa/.fasta for phmmer/jackhmmer/hmmsearch, .hmm for hmmsearch)
        --hmmer_tool          phmmer | jackhmmer | hmmsearch      [default: phmmer]
        --hmmer_evalue        Reporting E-value threshold (-E)    [default: 10.0]
        --hmmer_incE          Inclusion E-value threshold (--incE)[default: 0.001]
        --hmmer_N             jackhmmer max iterations (-N)       [default: 5]
        --hmmer_cpu           CPU threads (--cpu)                 [default: 2]
        --hmmer_use_original_alignment  Skip hmmalign re-alignment when building HMM from FASTA [default: false]

    Step modes (re-run a specific phase; prior steps must be complete):
        --extract             Re-run region extraction + ORF finding + allele_info (steps 1-5)

    Alternative step syntax: --step extract

    Profiles:
        -profile docker       Run with Docker (recommended)
        -profile standard     Run with local conda env

    Examples:
        # Full pipeline
        nextflow run main.nf \\
            --allele_search_dir ./allele_results/blast_summary_filtered \\
            --genome_dir /data/genomes \\
            --region_start 0 --region_end 500 --strand_mode both \\
            --outdir ./results -profile docker -resume

        # Re-run ORF finding + allele_info only
        nextflow run main.nf --extract \\
            --allele_search_dir ./allele_results/blast_summary_filtered \\
            --genome_dir /data/genomes \\
            --region_start 0 --region_end 500 \\
            --outdir ./results -profile docker -resume
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ---------------------------------------------------------------------------
// Step mode resolver
// ---------------------------------------------------------------------------
def stepMode = params.step ?: (
    params.extract ? 'extract' : 'full'
)

// ---------------------------------------------------------------------------
// Parameter validation
// ---------------------------------------------------------------------------
if (stepMode in ['full', 'extract']) {
    if (!params.allele_search_dir) error "--allele_search_dir is required"
    if (!params.genome_dir)        error "--genome_dir is required"
    if (params.region_start == null) error "--region_start is required"
    if (params.region_end   == null) error "--region_end is required"
}

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------
workflow {

    if (stepMode in ['full', 'extract']) {

        // Step 1: parse allele_search CSVs → hits.tsv
        PARSE_ALLELE_INPUT(Channel.value(file(params.allele_search_dir)))

        ch_genome_map = Channel.fromPath("${params.genome_dir}/*.{fasta,fa,fna,ffn}")
            .map { f -> [f.baseName, f] }

        // Step 2: group all hits for the same genome into one task (genome loaded once)
        ch_hits = PARSE_ALLELE_INPUT.out
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                tuple(row.species_name,
                      tuple(row.contig,
                            row.subject_start as int,
                            row.subject_end as int,
                            row.source_gene))
            }
            .groupTuple()                        // → [species_name, [[c,s,e,gene], ...]]
            .combine(ch_genome_map, by: 0)       // → [species_name, [[hits...]], genome_fasta]
            .map { sp, hits, genome -> tuple(sp, genome, hits) }

        // Batch N genomes per task to reduce task overhead (~7525 → ~75 tasks at default)
        ch_extract_batches = ch_hits
            .buffer(size: params.extract_batch_size, remainder: true)
            .map { batch ->
                tuple(
                    batch.collect { it[1] },   // list of genome Path objects → staged as path
                    batch.collect { it[0] },   // list of species name strings
                    batch.collect { it[2] }    // list of hits lists
                )
            }

        EXTRACT_DOWNSTREAM_REGION(ch_extract_batches)

        // Step 3: flatten per-genome multi-file output → batch 100 regions per FIND_ORFS task
        FIND_ORFS(
            EXTRACT_DOWNSTREAM_REGION.out
                .flatten()
                .buffer(size: params.batch_size, remainder: true)
        )

        // Step 4: collect all batched ORF FASTAs into one file
        COLLECT_ORFS(FIND_ORFS.out.collect())

        // Step 3.5 (optional): single HMMER run on unique proteins from all collected ORFs
        def hmmer_hits_ch
        if (params.hmmer_query) {
            HMMER_SEARCH(COLLECT_ORFS.out, file(params.hmmer_query))
            hmmer_hits_ch = HMMER_SEARCH.out[0]
        } else {
            hmmer_hits_ch = Channel.value([])
        }

        // Step 5: translate + group by protein → allele_info + per_genome_orfs output
        ORF_ALLELE_INFO(COLLECT_ORFS.out, hmmer_hits_ch)

    }
}

// ---------------------------------------------------------------------------
// Completion handler
// ---------------------------------------------------------------------------
workflow.onComplete {
    def missing = []

    if (!file("${params.outdir}/candidate_orfs").exists())
        missing << "candidate_orfs/  (COLLECT_ORFS produced no output)"
    if (!file("${params.outdir}/allele_info").exists())
        missing << "allele_info/  (ORF_ALLELE_INFO produced no output)"
    if (!file("${params.outdir}/per_genome_orfs").exists())
        missing << "per_genome_orfs/  (ORF_ALLELE_INFO produced no output)"

    if (missing) {
        log.warn "======================================================"
        log.warn " Pipeline finished but expected outputs are MISSING:"
        missing.each { log.warn "   - $it" }
        log.warn " Check .nextflow.log and the work/ directory for details."
        log.warn "======================================================"
    }
}
