#!/usr/bin/env nextflow

process CallVariants {

    container 'docker://google/deepvariant:1.10.0'

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(reference), path(fai)
    path regions_file

    output:
    tuple val(sample), path("${sample}.variants.vcf.gz"), emit: vcf
    // DeepVariant writes a standalone QC HTML (default --vcf_stats_report). MultiQC has no
    // DeepVariant module, so it is published as-is (not fed to MultiQC). optional: a user can
    // disable it via deepvariant_args (--vcf_stats_report=false).
    tuple val(sample), path("${sample}.variants.visual_report.html"), optional: true, emit: report

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    mkdir -p tmp

    # Read regions from file (space-separated chromosome names), trimmed
    REGIONS=\$(cat ${regions_file} | tr -s ' ' | sed 's/^ //;s/ \$//')

    # Build --regions only if regions were detected. Use a bash array so the whole
    # space-separated list is passed as ONE quoted value ("1 2 3 ... MT"); passing it
    # unquoted word-splits it and DeepVariant only keeps the first chromosome.
    REGIONS_ARGS=()
    if [ -n "\$REGIONS" ]; then
        REGIONS_ARGS+=(--regions "\$REGIONS")
    fi

    /opt/deepvariant/bin/run_deepvariant \\
        --model_type PACBIO \\
        --ref ${reference} \\
        --reads ${bam} \\
        --output_vcf ${sample}.variants.vcf.gz \\
        --intermediate_results_dir tmp \\
        --num_shards ${task.cpus} \\
        \${REGIONS_ARGS[@]+"\${REGIONS_ARGS[@]}"} \\
        ${args}
    """
}
