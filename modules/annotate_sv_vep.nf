#!/usr/bin/env nextflow

process AnnotateSVWithVEP {
    container 'ensemblorg/ensembl-vep:release_115.2'

    input:
    tuple val(sample), val(chromosome), val(interval), path(vcf), path(tbi), path(annotation), path(reference), path(fai)

    output:
    tuple val(sample), val(chromosome), val(interval), path("${sample}_sv_annotated_chr${chromosome}_${interval}.vcf.gz"), path("${sample}_sv_annotated_chr${chromosome}_${interval}.vcf.gz.tbi"), emit: annotated_sv
    tuple val(sample), val(chromosome), val(interval), path("${sample}_sv_summary_chr${chromosome}_${interval}.txt"), emit: summary
    path "*_summary.html", optional: true, emit: vep_summary

    script:
    def args      = task.ext.args ?: ''
    def use_gff   = params.vep_gff != null
    def annotated = "${sample}_sv_annotated_chr${chromosome}_${interval}.vcf"
    def summary   = "${sample}_sv_summary_chr${chromosome}_${interval}.txt"

    // GFF prep block - only emitted in GFF mode
    // Each VEP process prepares its own copy in its work dir (independent tasks).
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

    // Cache-only features (require pre-computed data not available via GFF). --sift is
    // dropped entirely here: it scores amino-acid substitutions, meaningless for SVs.
    def vep_cache_features = use_gff \
        ? "" \
        : "--gene_phenotype --regulatory"

    """
    set -euo pipefail

    # Number of SVs in the sliced interval VCF (skip header).
    N_SV=\$(zcat ${vcf} | grep -vc '^#' || true)
    echo "Number of SVs in input: \$N_SV"

    if [ "\$N_SV" -eq 0 ]; then
        echo "WARNING: No SVs in interval chr${chromosome}:${interval} - emitting header-only annotated VCF"
        zcat ${vcf} > ${annotated}
    else
${gff_prep}
        # VEP on structural variants: it reads symbolic alleles (<DEL>/<DUP>/<INS>/<INV>)
        # and BND records, using POS + INFO/END + INFO/SVLEN to define the affected span,
        # then overlaps the SAME gene model as the SNV path and assigns SV consequence terms
        # (transcript_ablation, feature_truncation, ...). SNV-only scoring (--sift) is dropped.
        vep \\
            --input_file ${vcf} \\
            --output_file ${annotated} \\
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

    # Build the per-interval SV summary FROM the annotated VCF so it carries VEP gene/
    # consequence. CSQ order (VEP --symbol): Allele|Consequence|IMPACT|SYMBOL|... so the
    # consequence is field 2 and the gene SYMBOL is field 4.
    awk -F'\\t' '
    /^#/ { next }
    {
        nsv++
        csq=""; svtype="."; svlen="."; endp="."
        n=split(\$8, kv, ";")
        for (j=1; j<=n; j++) {
            split(kv[j], p, "=")
            if (p[1]=="SVTYPE") svtype=p[2]
            else if (p[1]=="SVLEN") svlen=p[2]
            else if (p[1]=="END") endp=p[2]
            else if (p[1]=="CSQ") csq=p[2]
        }
        gt="."
        if (NF>=10) { split(\$10, g, ":"); gt=g[1] }
        gene="."; cons="."
        if (csq != "") {
            split(csq, tr, ",")
            split(tr[1], f, "|")
            cons = (f[2]=="" ? "." : f[2])
            gene = (f[4]=="" ? "." : f[4])
        }
        types[svtype]++
        row[nsv] = sprintf("%-6s  %-12s  %-12s  %-9s  %-8s  %-6s  %-14s  %-30s", svtype, \$2, endp, svlen, \$7, gt, gene, cons)
    }
    END {
        if (nsv == 0) {
            print "No structural variants found overlapping this region."
        } else {
            printf "Structural variants overlapping region: %d\\n\\n", nsv
            printf "By type:\\n"
            for (t in types) printf "  - %-6s %d\\n", t, types[t]
            printf "\\n"
            printf "%-6s  %-12s  %-12s  %-9s  %-8s  %-6s  %-14s  %-30s\\n", "TYPE", "POS", "END", "SVLEN", "FILTER", "GT", "GENE", "CONSEQUENCE"
            printf "%-6s  %-12s  %-12s  %-9s  %-8s  %-6s  %-14s  %-30s\\n", "------", "------------", "------------", "---------", "--------", "------", "--------------", "------------------------------"
            for (i=1; i<=nsv; i++) print row[i]
            printf "\\nNote: GENE/CONSEQUENCE = VEP annotation from the first transcript only.\\n"
        }
    }' ${annotated} > ${summary}

    # Compress + index the annotated SV VCF (both branches produced a plain .vcf).
    bgzip ${annotated}
    tabix -p vcf ${annotated}.gz
    """
}
