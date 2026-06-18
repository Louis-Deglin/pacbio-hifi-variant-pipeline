#!/usr/bin/env nextflow

process Align {

    container 'community.wave.seqera.io/library/pbmm2:26.1.99--118eceb01c386181'

    input:
    tuple val(sample), path(input)
    path reference

    output:
    tuple val(sample), path("${sample}.aligned.bam"), path("${sample}.aligned.bam.bai"), emit: bam

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    pbmm2 align \\
        ${reference} \\
        ${input} \\
        ${sample}.aligned.bam \\
        -j ${task.cpus} \\
        --sort \\
        ${args}
    """
}
