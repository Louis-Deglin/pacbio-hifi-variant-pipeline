#!/usr/bin/env nextflow

process CallSV {
    container 'community.wave.seqera.io/library/sawfish:2.2.1--884a5d05a2b2005f'

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(reference), path(fai)

    output:
    tuple val(sample), path("${sample}.sv.vcf.gz"), emit: sv_vcf

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    # sawfish is a two-step PacBio HiFi structural-variant caller:
    #   1) discover    -> per-sample candidate SVs + local assembly evidence (writes a dir)
    #   2) joint-call  -> genotyped VCF from one or more discover dirs
    sawfish discover \\
        --threads ${task.cpus} \\
        --ref ${reference} \\
        --bam ${bam} \\
        --output-dir discover_${sample} \\
        ${args}

    sawfish joint-call \\
        --threads ${task.cpus} \\
        --sample discover_${sample} \\
        --output-dir joint_${sample}

    # sawfish writes joint_<sample>/genotyped.sv.vcf.gz (VCF 4.4). Rename to a
    # sample-prefixed name to avoid collisions across parallel tasks. SliceSV
    # (bcftools) indexes it; no tabix is assumed inside the sawfish container.
    cp joint_${sample}/genotyped.sv.vcf.gz ${sample}.sv.vcf.gz
    """
}
