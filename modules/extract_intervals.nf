#!/usr/bin/env nextflow

process ExtractIntervals {
    container 'community.wave.seqera.io/library/bcftools:1.23.1--4d193a5f61d4aed7'

    input:
    tuple val(sample), path(vcf), path(tbi), val(chromosome), val(interval)

    output:
    tuple val(sample), val(chromosome), val(interval), path("${sample}_interval_chr${chromosome}_${interval}.vcf.gz"), path("${sample}_interval_chr${chromosome}_${interval}.vcf.gz.tbi"), emit: interval_vcf

    script:
    """
    set -euo pipefail

    # Input VCF (phased) already comes with its .tbi index from HiPhase.

    # Extract variants for this interval
    bcftools view -r ${chromosome}:${interval} ${vcf} -O z -o ${sample}_interval_chr${chromosome}_${interval}.vcf.gz

    # Index the interval VCF
    bcftools index -t ${sample}_interval_chr${chromosome}_${interval}.vcf.gz
    """
}
