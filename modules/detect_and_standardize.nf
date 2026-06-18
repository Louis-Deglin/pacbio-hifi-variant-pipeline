#!/usr/bin/env nextflow

process DetectAndStandardize {

    container 'community.wave.seqera.io/library/samtools:1.23.1--d76a06ff3aefee52'

    input:
    path reference

    output:
    tuple path("reference_std.fa"), path("reference_std.fa.fai"), emit: standardized
    path("regions.txt"), emit: regions

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # 1. Detect format by looking at the first sequence header
    FIRST_HEADER=\$(grep -m1 '^>' ${reference})
    FIRST_SEQ=\$(echo "\$FIRST_HEADER" | sed 's/>//' | awk '{print \$1}')

    if [[ "\$FIRST_SEQ" == NC_* ]]; then
        FORMAT="ncbi"
        echo "Detected NCBI format (first sequence: \$FIRST_SEQ)"
    else
        FORMAT="ensembl"
        echo "Detected Ensembl format (first sequence: \$FIRST_SEQ)"
    fi

    # 2. Standardize the FASTA
    if [[ "\$FORMAT" == "ncbi" ]]; then
        echo "Renaming NCBI headers to simple chromosome names..."

        awk '
        /^>/ {
            # Skip scaffolds (NW_, NT_, NZ_)
            if (\$0 ~ /^>(NW_|NT_|NZ_)/) {
                skip = 1
                next
            }
            skip = 0

            # Extract chromosome name using POSIX-compliant match()
            if (match(\$0, /chromosome [0-9A-Za-z]+,/)) {
                chrom = substr(\$0, RSTART, RLENGTH)
                sub(/^chromosome /, "", chrom)
                sub(/,\$/, "", chrom)
                print ">" chrom
            } else if (\$0 ~ /mitochondrion/) {
                print ">MT"
            } else {
                # Unknown NC_ sequence - keep original name as fallback
                split(\$0, parts, " ")
                print parts[1]
            }
            next
        }
        {
            if (skip == 0) print
        }
        ' ${reference} > reference_std.fa

        echo "FASTA standardized (NCBI -> simple names)"

    else
        # Ensembl: create symlink (no modification needed)
        ln -s \$(readlink -f ${reference}) reference_std.fa
        echo "Ensembl format: symlink created"
    fi

    # 3. Index the standardized FASTA
    samtools faidx reference_std.fa

    # 4. Extract chromosome regions from .fai
    grep -E '^([0-9]+|X|Y|Z|W|MT)\\s' reference_std.fa.fai | cut -f1 | tr '\\n' ' ' > regions.txt

    echo "Detected regions: \$(cat regions.txt)"
    """
}
