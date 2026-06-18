#!/usr/bin/env nextflow

process MultiQC {
    container 'community.wave.seqera.io/library/multiqc:1.35--5f40ae3381b1c04b'

    input:
    path qc_files

    output:
    tuple path("multiqc_report.html"), path("multiqc_data"), emit: report

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    multiqc . ${args}
    """
}
