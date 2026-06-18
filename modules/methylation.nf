#!/usr/bin/env nextflow

process MethylationProfiling {
    container 'community.wave.seqera.io/library/pb-cpg-tools:3.0.0--29244e31efc390ff'

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(reference), path(fai)

    output:
    tuple val(sample), path("${sample}.*.bed.gz"), optional: true, emit: beds
    tuple val(sample), path("${sample}.*.bw"), optional: true, emit: bigwig
    tuple val(sample), path("${sample}.no_methylation.txt"), optional: true, emit: skipped

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    # pb-CpG-tools reads 5mC from the MM/ML tags carried natively by PacBio HiFi BAMs and
    # writes per-CpG methylation scores. With a HiPhase haplotagged BAM (HP tags) it also
    # emits hap1/hap2 tracks (allele-specific methylation). Output naming with
    # --output-prefix <sample>: <sample>.combined.bed.gz (+ hap1/hap2) and matching .bw.
    #
    # Wrapped in `if ... then` so a BAM lacking MM/ML tags does NOT abort the process
    # (set -e would otherwise kill the run on a non-zero exit).
    if aligned_bam_to_cpg_scores \\
        --bam ${bam} \\
        --ref ${reference} \\
        --output-prefix ${sample} \\
        --modsites-mode reference \\
        --threads ${task.cpus} \\
        ${args}
    then
        :
    fi

    # If no methylation BED was produced (MM/ML tags absent), emit a sentinel so the run
    # continues and the per-interval report shows "No methylation profile found.".
    if [ ! -f ${sample}.combined.bed.gz ]; then
        echo "No methylation profile found." > ${sample}.no_methylation.txt
    fi
    """
}
