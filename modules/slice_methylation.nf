#!/usr/bin/env nextflow

process SliceMethylation {
    container 'community.wave.seqera.io/library/bcftools:1.23.1--4d193a5f61d4aed7'

    input:
    tuple val(sample), path(beds), val(chromosome), val(interval)

    output:
    tuple val(sample), val(chromosome), val(interval), path("${sample}_*_chr${chromosome}_${interval}.bed"), optional: true, emit: sliced
    tuple val(sample), val(chromosome), val(interval), path("${sample}_methylation_summary_chr${chromosome}_${interval}.txt"), emit: summary

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # pb-CpG-tools BED columns: 1=chrom 2=start(0-based) 3=end 4=meth_percent 5=haplotype 6=coverage
    # Slice each available track to the region with bgzip+awk (no tabix dependency). A CpG
    # [start,end) overlaps the 1-based region s-e when start<=e and end>=s.
    for track in combined hap1 hap2; do
        bed="${sample}.\${track}.bed.gz"
        out="${sample}_\${track}_chr${chromosome}_${interval}.bed"
        if [ -f "\$bed" ]; then
            bgzip -dc "\$bed" \\
                | awk -v c="${chromosome}" -v iv="${interval}" \\
                    'BEGIN{split(iv,r,"-"); s=r[1]; e=r[2]} \$1==c && \$2<=e && \$3>=s' \\
                > "\$out" || true
        fi
    done

    COMBINED="${sample}_combined_chr${chromosome}_${interval}.bed"
    HAP1="${sample}_hap1_chr${chromosome}_${interval}.bed"
    HAP2="${sample}_hap2_chr${chromosome}_${interval}.bed"
    SUMMARY="${sample}_methylation_summary_chr${chromosome}_${interval}.txt"

    # Ensure hap files exist so the awk getline below never trips on a missing file.
    [ -f "\$HAP1" ] || touch "\$HAP1"
    [ -f "\$HAP2" ] || touch "\$HAP2"

    if [ ! -s "\$COMBINED" ]; then
        echo "No methylation profile found." > "\$SUMMARY"
    else
        # Summary stats (CpG count, mean % combined/hap1/hap2) + per-CpG table.
        # hap1/hap2 percentages are matched to combined sites by start coordinate.
        awk -v hap1f="\$HAP1" -v hap2f="\$HAP2" '
        BEGIN {
            while ((getline line < hap1f) > 0) { n=split(line,a,"\\t"); h1[a[2]]=a[4] }
            while ((getline line < hap2f) > 0) { n=split(line,a,"\\t"); h2[a[2]]=a[4] }
        }
        {
            start[NR]=\$2; pos[NR]=\$3; meth[NR]=\$4; cov[NR]=\$6
            sumc+=\$4; nc++
        }
        END {
            mc = (nc>0) ? sumc/nc : 0
            s1=0; n1=0; for (k in h1){ s1+=h1[k]; n1++ }
            s2=0; n2=0; for (k in h2){ s2+=h2[k]; n2++ }

            printf "CpG sites in region:   %d\\n", nc
            printf "Mean methylation:      %.1f %%   (combined)\\n", mc
            if (n1>0) printf "  - Haplotype 1:        %.1f %%\\n", s1/n1; else printf "  - Haplotype 1:        N/A\\n"
            if (n2>0) printf "  - Haplotype 2:        %.1f %%\\n", s2/n2; else printf "  - Haplotype 2:        N/A\\n"
            printf "\\n"
            printf "%-10s  %-7s  %-5s  %-7s  %-7s\\n", "POS", "METH%", "COV", "HAP1%", "HAP2%"
            printf "%-10s  %-7s  %-5s  %-7s  %-7s\\n", "----------", "-------", "-----", "-------", "-------"
            for (i=1; i<=NR; i++) {
                v1 = (start[i] in h1) ? sprintf("%.1f", h1[start[i]]) : "."
                v2 = (start[i] in h2) ? sprintf("%.1f", h2[start[i]]) : "."
                printf "%-10s  %-7.1f  %-5s  %-7s  %-7s\\n", pos[i], meth[i], cov[i], v1, v2
            }
        }' "\$COMBINED" > "\$SUMMARY"
    fi
    """
}
