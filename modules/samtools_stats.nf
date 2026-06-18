#!/usr/bin/env nextflow

process SamtoolsStats {
    container 'community.wave.seqera.io/library/samtools:1.23.1--d76a06ff3aefee52'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path "${sample}.samtools_stats.txt", emit: stats

    script:
    """
    set -euo pipefail

    samtools stats -@ ${task.cpus} ${bam} > ${sample}.samtools_stats.txt
    """
}
