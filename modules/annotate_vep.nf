#!/usr/bin/env nextflow

process AnnotateWithVEP {
    container 'ensemblorg/ensembl-vep:release_115.2'

    input:
    tuple val(sample), val(chromosome), val(interval), path(vcf), path(tbi), path(annotation), path(reference), path(fai)

    output:
    tuple val(sample), val(chromosome), val(interval), path("${sample}_interval_chr${chromosome}_${interval}_annotated.vcf.gz"), path("${sample}_interval_chr${chromosome}_${interval}_annotated.vcf.gz.tbi"), emit: annotated_vcf
    path "*_summary.html", optional: true, emit: vep_summary

    script:
    def args     = task.ext.args ?: ''
    def use_gff  = params.vep_gff != null

    // GFF prep block - only emitted in GFF mode. Follows the official Ensembl recipe
    // (https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html):
    // strip comments, sort by chrom + start + end with explicit TAB separator
    // (-t \$'\\t' is critical, otherwise attribute fields with spaces break the sort),
    // then bgzip and tabix.
    def gff_prep = use_gff ? """
        zcat -f ${annotation} \\
            | grep -v "^#" \\
            | sort -k1,1 -k4,4n -k5,5n -t \$'\\t' \\
            | bgzip -c > prepared_annotation.gff.gz
        tabix -p gff prepared_annotation.gff.gz
""" : ""

    // Annotation source flags (cache vs gff)
    def vep_source = use_gff \
        ? "--gff prepared_annotation.gff.gz --fasta ${reference}" \
        : "--cache --dir_cache ${annotation} --fasta ${reference} --offline"

    // Cache-only features (require pre-computed data not available via GFF)
    def vep_cache_features = use_gff \
        ? "" \
        : "--sift b --gene_phenotype --regulatory"

    """
    set -euo pipefail

    # Check if VCF contains any variants (skip header lines)
    N_VARIANTS=\$(zcat ${vcf} | grep -vc '^#' || true)
    echo "Number of variants in input: \$N_VARIANTS"

    if [ "\$N_VARIANTS" -eq 0 ]; then
        echo "WARNING: No variants in interval chr${chromosome}:${interval} - emitting empty annotated VCF"
        # Copy the header-only VCF as the "annotated" output, then compress & index
        zcat ${vcf} > ${sample}_interval_chr${chromosome}_${interval}_annotated.vcf
    else
${gff_prep}
        # Run VEP annotation
        vep \\
            --input_file ${vcf} \\
            --output_file ${sample}_interval_chr${chromosome}_${interval}_annotated.vcf \\
            --vcf \\
            --species ${params.vep_species} \\
            --assembly ${params.vep_assembly} \\
            ${vep_source} \\
            --fork ${task.cpus} \\
            --symbol \\
            --canonical \\
            --biotype \\
            ${vep_cache_features} \\
            --force_overwrite \\
            ${args}
    fi

    # Compress and index annotated VCF (both branches produce a plain .vcf)
    bgzip ${sample}_interval_chr${chromosome}_${interval}_annotated.vcf
    tabix -p vcf ${sample}_interval_chr${chromosome}_${interval}_annotated.vcf.gz
    """
}
