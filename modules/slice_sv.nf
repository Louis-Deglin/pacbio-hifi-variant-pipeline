#!/usr/bin/env nextflow

process SliceSV {
    container 'community.wave.seqera.io/library/bcftools:1.23.1--4d193a5f61d4aed7'

    input:
    tuple val(sample), path(sv_vcf), val(chromosome), val(interval)

    output:
    tuple val(sample), val(chromosome), val(interval), path("${sample}_sv_chr${chromosome}_${interval}.vcf.gz"), path("${sample}_sv_chr${chromosome}_${interval}.vcf.gz.tbi"), emit: sliced

    script:
    """
    #!/bin/bash
    set -euo pipefail

    OUT="${sample}_sv_chr${chromosome}_${interval}.vcf.gz"

    # Region queries need an index on the genome-wide sawfish VCF.
    bcftools index -t -f ${sv_vcf}

    # --regions-overlap 2: keep any SV that OVERLAPS the region. An SV (DEL/DUP/INV...)
    # can extend beyond the interval without being fully contained, unlike a point
    # variant, so containment-only slicing would miss breakpoints sitting on the edge.
    # The sliced VCF is then annotated by AnnotateSVWithVEP 
    bcftools view --regions-overlap 2 -r ${chromosome}:${interval} ${sv_vcf} -Oz -o "\$OUT"
    bcftools index -t -f "\$OUT"
    """
}
