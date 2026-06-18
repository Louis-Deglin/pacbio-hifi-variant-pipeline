#!/usr/bin/env nextflow

process HiPhase {
    container 'community.wave.seqera.io/library/bcftools_hiphase_samtools:7daa1f30f0e0c390'

    input:
    tuple val(sample), path(bam), path(bai), path(vcf)
    tuple path(reference), path(fai)

    output:
    tuple val(sample), path("${sample}.phased.vcf.gz"), path("${sample}.phased.vcf.gz.tbi"), emit: phased_vcf
    tuple val(sample), path("${sample}.haplotagged.bam"), path("${sample}.haplotagged.bam.bai"), emit: haplotagged_bam
    tuple val(sample), path("${sample}.phase_stats.tsv"), path("${sample}.phase_summary.tsv"), path("${sample}.phase_blocks.tsv"), emit: stats

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    # HiPhase requires an indexed input VCF (CallVariants emits only the .vcf.gz).
    bcftools index -t ${vcf}

    # Phase variants and haplotag the BAM. --ignore-read-groups avoids matching the
    # BAM RG/SM tag to the VCF sample name (single sample per BAM here).
    hiphase \\
        --bam ${bam} \\
        --vcf ${vcf} \\
        --output-vcf ${sample}.phased.vcf.gz \\
        --output-bam ${sample}.haplotagged.bam \\
        --reference ${reference} \\
        --threads ${task.cpus} \\
        --ignore-read-groups \\
        --stats-file ${sample}.phase_stats.tsv \\
        --summary-file ${sample}.phase_summary.tsv \\
        --blocks-file ${sample}.phase_blocks.tsv \\
        ${args}

    # HiPhase already indexes its own outputs (${sample}.phased.vcf.gz.tbi and
    # ${sample}.haplotagged.bam.bai) - re-indexing here would fail (.tbi exists).
    """
}
