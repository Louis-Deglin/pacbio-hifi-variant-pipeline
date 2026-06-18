#!/usr/bin/env nextflow

// Run-level aggregated HTML, built ONCE from every per-interval fragment:
//   - general_summary.html : all (sample x interval) fragments + a plain-text table
//     of contents (the aggregated report covering every sample and interval).
//
// Reuses the MultiQC container (python3). The build script and the manifest are
// staged as inputs; the fragments are collected from IntervalReport.

process RunReport {
    container 'community.wave.seqera.io/library/multiqc:1.35--5f40ae3381b1c04b'

    input:
    path fragments
    path manifest
    path build_script

    output:
    path "general_summary.html", emit: summary

    script:
    def run_name = params.run_name ?: ''
    """
    set -euo pipefail

    python3 ${build_script} \\
        --manifest '${manifest}' \\
        --run-name '${run_name}' \\
        --out-summary general_summary.html
    """
}
