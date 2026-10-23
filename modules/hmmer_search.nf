// SHP — hmmer_search.nf
// Process: HMMER_SEARCH
//
// Step 3.5 (optional): Single HMMER run on all unique proteins from the full collected
//   ORF set. Skipped if --hmmer_query is not provided.
//
// Runs after COLLECT_ORFS (not per-region). Translates all ORFs, deduplicates by protein
// sequence, runs HMMER once on unique proteins, then expands the tblout so every ORF
// occurrence of a hit protein appears — preserving full_header matching in orf_allele_info.py.
//
// Query file detection:
//   .hmm extension  → hmmsearch with the pre-built profile directly
//   FASTA extension → phmmer / jackhmmer / hmmsearch (depending on --hmmer_tool)
//                     For hmmsearch: hmmbuild [+ hmmalign + hmmbuild] to build profile first

process HMMER_SEARCH {
    publishDir "${params.outdir}/hmmer", mode: 'copy'
    stageInMode 'copy'

    input:
    path raw_orfs    // all_orfs_raw.fasta from COLLECT_ORFS
    path query_file

    output:
    path "all_hits.tblout"
    path "jackhmmer_final.hmm", optional: true

    script:
    // ---------- Groovy: resolve tool and build command fragments ----------
    def is_hmm      = query_file.extension == 'hmm'
    def common_opts = "-E ${params.hmmer_evalue} --incE ${params.hmmer_incE} --cpu ${params.hmmer_cpu} --tblout raw_hits.tblout --noali"

    def pre_search_cmds
    def search_cmd

    if (is_hmm) {
        pre_search_cmds = ""
        search_cmd      = "hmmsearch ${common_opts} ${query_file} unique_proteins.faa"

    } else if (params.hmmer_tool == "phmmer") {
        pre_search_cmds = ""
        search_cmd      = "phmmer ${common_opts} ${query_file} unique_proteins.faa"

    } else if (params.hmmer_tool == "jackhmmer") {
        pre_search_cmds = ""
        search_cmd      = "jackhmmer -N ${params.hmmer_N} --chkhmm jackhmmer_ckpt ${common_opts} ${query_file} unique_proteins.faa"

    } else {
        // hmmsearch with FASTA query: build HMM profile, optionally refine with hmmalign
        if (params.hmmer_use_original_alignment) {
            pre_search_cmds = "hmmbuild query.hmm ${query_file}"
        } else {
            pre_search_cmds = """\
hmmbuild initial.hmm ${query_file}
hmmalign --trim --outformat afa initial.hmm ${query_file} > realigned.afa
hmmbuild query.hmm realigned.afa"""
        }
        search_cmd = "hmmsearch ${common_opts} query.hmm unique_proteins.faa"
    }
    // ---------- end Groovy ----------

    """
    # Translate all ORFs; deduplicate by protein sequence; write unique_proteins.faa.
    # Build protein → [all full_headers] map so tblout can be expanded below.
    python3 - "${raw_orfs}" <<'PYEOF'
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

protein_to_headers = defaultdict(list)
with open('unique_proteins.faa', 'w') as out:
    for rec in SeqIO.parse(sys.argv[1], 'fasta'):
        dna = str(rec.seq).upper().replace('-', '').replace('N', '')
        dna = dna[:len(dna) - len(dna) % 3]
        if len(dna) < 3:
            continue
        prot = str(Seq(dna).translate())
        if prot.endswith('*'):
            prot = prot[:-1]
        if not prot:
            continue
        prot = 'M' + prot[1:]
        protein_to_headers[prot].append(rec.description)
        if len(protein_to_headers[prot]) == 1:   # first occurrence → write representative
            out.write('>' + rec.description + '\\n' + prot + '\\n')
PYEOF

    # Run HMMER once on unique proteins (guard: skip if nothing translatable)
    if [ -s unique_proteins.faa ]; then
        ${pre_search_cmds}
        ${search_cmd}
    else
        touch raw_hits.tblout
    fi

    # Save final jackhmmer HMM checkpoint (only present when --hmmer_tool jackhmmer)
    last_ckpt=\$(ls jackhmmer_ckpt-*.hmm 2>/dev/null | sort -t- -k2 -n | tail -1)
    if [ -n "\$last_ckpt" ]; then
        cp "\$last_ckpt" jackhmmer_final.hmm
    fi

    # Expand tblout: for each hit representative, emit one line per ORF occurrence
    # of that protein so orf_allele_info.py full_header matching works unchanged.
    python3 - "${raw_orfs}" <<'PYEOF2'
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

# Rebuild protein → [all headers] map (same translation logic as above)
protein_to_headers = defaultdict(list)
for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    dna = str(rec.seq).upper().replace('-', '').replace('N', '')
    dna = dna[:len(dna) - len(dna) % 3]
    if len(dna) < 3:
        continue
    prot = str(Seq(dna).translate())
    if prot.endswith('*'):
        prot = prot[:-1]
    if not prot:
        continue
    prot = 'M' + prot[1:]
    protein_to_headers[prot].append(rec.description)

# Map representative header (first occurrence) → protein sequence
rep_to_protein = {hdrs[0]: prot for prot, hdrs in protein_to_headers.items()}

with open('raw_hits.tblout') as raw, open('all_hits.tblout', 'w') as out:
    for line in raw:
        if line.startswith('#') or not line.strip():
            out.write(line)
            continue
        parts = line.split()
        rep_header = parts[0]
        if rep_header in rep_to_protein:
            prot = rep_to_protein[rep_header]
            for full_hdr in protein_to_headers[prot]:
                parts[0] = full_hdr
                out.write('\\t'.join(parts) + '\\n')
        else:
            out.write(line)
PYEOF2
    """
}
