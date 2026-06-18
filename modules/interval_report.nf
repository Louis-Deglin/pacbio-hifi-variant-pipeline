#!/usr/bin/env nextflow

// Per-(sample x interval) HTML report.
// Parses the annotated VCF and (when --run_sv) the annotated SV VCF, embeds the
// methylation summary, and renders a rich HTML page (curated CSQ columns + a
// collapsible full-CSQ <details> per variant). Also emits a reusable HTML fragment
// that RunReport concatenates into the run-level full_report.html.
//
// Reuses the MultiQC container (it ships python3); the render script is staged as an
// input (file from projectDir) so no executable bit / PATH magic is needed.

process IntervalReport {
    container 'community.wave.seqera.io/library/multiqc:1.35--5f40ae3381b1c04b'

    input:
    tuple val(sample), val(chromosome), val(interval), path(snv_vcf), path(meth_summary), path(sv_vcf)
    path render_script

    output:
    tuple val(sample), val(chromosome), val(interval), path("${sample}_report_chr${chromosome}_${interval}.html"), emit: report
    path "${sample}_fragment_chr${chromosome}_${interval}.html", emit: fragment

    script:
    def run_meth = params.run_methylation ? 'true' : 'false'
    def run_sv   = params.run_sv ? 'true' : 'false'
    def standalone = "${sample}_report_chr${chromosome}_${interval}.html"
    def fragment   = "${sample}_fragment_chr${chromosome}_${interval}.html"
    """
    set -euo pipefail

    python3 ${render_script} \\
        --sample '${sample}' \\
        --chrom '${chromosome}' \\
        --interval '${interval}' \\
        --snv-vcf '${snv_vcf}' \\
        --sv-vcf '${sv_vcf}' \\
        --meth-summary '${meth_summary}' \\
        --run-sv ${run_sv} \\
        --run-methylation ${run_meth} \\
        --out-standalone '${standalone}' \\
        --out-fragment '${fragment}'
    """
}
