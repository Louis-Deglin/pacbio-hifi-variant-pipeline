#!/usr/bin/env nextflow

process BcftoolsStats {
    container 'community.wave.seqera.io/library/bcftools:1.23.1--4d193a5f61d4aed7'

    input:
    tuple val(sample), path(vcf), val(tag)

    output:
    path "*.bcftools_stats.txt", emit: stats

    script:
    // tag distinguishes call sets when the process is reused via an alias (e.g. 'sv').
    // Empty tag keeps the original <sample>.bcftools_stats.txt name.
    def suffix = tag ? ".${tag}" : ""
    """
    set -euo pipefail

    bcftools stats ${vcf} > ${sample}${suffix}.bcftools_stats.txt
    """
}
