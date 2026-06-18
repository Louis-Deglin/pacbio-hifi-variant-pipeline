#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =============================================================================
# Helper: for each tool, propose (1) show current flags, (2) prompt extra args.
# Result returned in $EXTRA_ARGS_RETURN.
# =============================================================================
EXTRA_ARGS_RETURN=""

ask_extra_args() {
    local tool_name="$1"
    local current_flags="$2"

    EXTRA_ARGS_RETURN=""

    echo ""
    echo "---- ${tool_name} ----"

    read -p "Show options currently passed to ${tool_name}? (y/n) [n]: " ANS
    if [ "${ANS:-n}" = "y" ]; then
        echo ""
        echo "Pipeline currently calls ${tool_name} with:"
        printf "  %s\n" "$current_flags"
        echo ""
    fi

    read -p "Extra ${tool_name} args (empty = none): " EXTRA_ARGS_RETURN
}

# =============================================================================
# Input mode: a single sample (interactive) OR several samples (CSV samplesheet)
# =============================================================================

echo "=== Input ==="
echo "  1  = a single sample  (provide a BAM + interval(s) interactively)"
echo "  2+ = multiple samples (provide a CSV samplesheet)"
echo ""
echo "CSV samplesheet format (columns): sample,bam,chromosome,interval"
echo "  - one row per (sample, interval); BAM paths must be ABSOLUTE."
echo "  - for several intervals of the SAME sample, repeat the row with the"
echo "    same sample + bam and a different chromosome/interval."
echo ""
read -p "How many samples? (1 / 2+) [1]: " SAMPLE_COUNT
SAMPLE_COUNT=${SAMPLE_COUNT:-1}

SAMPLESHEET=""        # set directly (multi mode) or generated (single mode)
SAMPLE_NAME=""
INPUT_BAM=""
CHROMS=()
INTERVALS=()
NUM_INTERVALS=0

if [ "$SAMPLE_COUNT" = "1" ]; then
    MODE="single"

    read -p "Sample name (ex: COW_001): " SAMPLE_NAME
    [ -z "$SAMPLE_NAME" ] && echo "ERROR: Sample name required" && exit 1

    read -p "Input BAM: " INPUT_BAM
    [ -z "$INPUT_BAM" ] && echo "ERROR: Input BAM required" && exit 1
else
    MODE="multi"

    read -p "Samplesheet CSV path: " SAMPLESHEET
    [ -z "$SAMPLESHEET" ] && echo "ERROR: CSV path required" && exit 1
    [ ! -f "$SAMPLESHEET" ] && echo "ERROR: CSV not found: $SAMPLESHEET" && exit 1
fi

read -p "Reference FASTA: " REFERENCE
[ -z "$REFERENCE" ] && echo "ERROR: Reference required" && exit 1

# =============================================================================
# VEP annotation source: cache directory OR GFF file (mutually exclusive)
# =============================================================================

echo ""
echo "=== VEP annotation source ==="
echo "  1) VEP cache directory (recommended when available, full annotations)"
echo "  2) GFF file"
echo ""
echo "NOTE: GFF mode disables the following VEP features (they require the cache):"
echo "        --sift           (SIFT pathogenicity prediction)"
echo "        --gene_phenotype (gene phenotype associations)"
echo "        --regulatory     (Ensembl regulatory features)"
echo "      Core consequence calling, --symbol, --canonical and --biotype still work"
echo "      provided the GFF declares them."
echo ""
read -p "Annotation source (1/2) [1]: " VEP_SOURCE_CHOICE
VEP_SOURCE_CHOICE=${VEP_SOURCE_CHOICE:-1}

VEP_CACHE=""
VEP_GFF=""

case "$VEP_SOURCE_CHOICE" in
    1)
        read -p "VEP cache directory: " VEP_CACHE
        [ -z "$VEP_CACHE" ] && echo "ERROR: VEP cache directory is required" && exit 1
        if [ ! -d "$VEP_CACHE" ]; then
            echo "ERROR: VEP cache directory not found: $VEP_CACHE"
            exit 1
        fi
        ;;
    2)
        read -p "GFF file (.gff or .gff.gz): " VEP_GFF
        [ -z "$VEP_GFF" ] && echo "ERROR: GFF file is required" && exit 1
        [ ! -f "$VEP_GFF" ] && echo "ERROR: GFF file not found: $VEP_GFF" && exit 1
        ;;
    *)
        echo "ERROR: Invalid choice (expected 1 or 2)"
        exit 1
        ;;
esac

# Sanity check: cache and GFF must be mutually exclusive, and one of them is required
if [ -n "$VEP_CACHE" ] && [ -n "$VEP_GFF" ]; then
    echo "ERROR: --vep_cache and --vep_gff are mutually exclusive; please provide only one"
    exit 1
fi
if [ -z "$VEP_CACHE" ] && [ -z "$VEP_GFF" ]; then
    echo "ERROR: you must provide either a VEP cache directory or a GFF file"
    exit 1
fi

# =============================================================================
# Species & assembly (passed to VEP)
# =============================================================================

echo ""
echo "=== Species & assembly ==="
read -p "VEP species [bos_taurus]: " VEP_SPECIES
VEP_SPECIES=${VEP_SPECIES:-bos_taurus}
read -p "VEP assembly [ARS-UCD2.0]: " VEP_ASSEMBLY
VEP_ASSEMBLY=${VEP_ASSEMBLY:-ARS-UCD2.0}

# =============================================================================
# Output directory
# =============================================================================

read -p "Output directory [Results]: " OUTPUT_DIR
OUTPUT_DIR=${OUTPUT_DIR:-Results}

# Run name: results land under <outdir>/<run_name>/ (one self-contained folder per run).
read -p "Run name: " RUN_NAME
[ -z "$RUN_NAME" ] && echo "ERROR: Run name required" && exit 1
RUN_DIR="${OUTPUT_DIR}/${RUN_NAME}"

# =============================================================================
# Intervals (single-sample mode only - multi mode reads them from the CSV)
# =============================================================================

if [ "$MODE" = "single" ]; then
    read -p "Multiple intervals? (y/n) [n]: " MULTIPLE
    MULTIPLE=${MULTIPLE:-n}

    if [ "$MULTIPLE" = "y" ]; then
        read -p "Number of intervals: " NUM_INTERVALS
        if ! [[ "$NUM_INTERVALS" =~ ^[0-9]+$ ]] || [ "$NUM_INTERVALS" -lt 1 ]; then
            echo "ERROR: Invalid number"
            exit 1
        fi
        for i in $(seq 1 "$NUM_INTERVALS"); do
            echo ""
            echo "=== Interval $i/$NUM_INTERVALS ==="
            read -p "Chromosome (ex: 26): " INTERVAL_CHROM
            read -p "Interval (ex: 14404993-14405000): " INTERVAL
            CHROMS+=("$INTERVAL_CHROM")
            INTERVALS+=("$INTERVAL")
        done
    else
        read -p "Chromosome (ex: 26): " INTERVAL_CHROM
        read -p "Interval (ex: 14404993-14405000): " INTERVAL
        CHROMS=("$INTERVAL_CHROM")
        INTERVALS=("$INTERVAL")
        NUM_INTERVALS=1
    fi
fi

# =============================================================================
# Optional extra arguments for each tool
# =============================================================================

echo ""
echo "=== Tool configuration (optional) ==="
read -p "Configure extra arguments for any tool? (y/n) [n]: " CONFIGURE_EXTRA
CONFIGURE_EXTRA=${CONFIGURE_EXTRA:-n}

ALIGN_ARGS=""
DEEPVARIANT_ARGS=""
VEP_ARGS=""

if [ "$CONFIGURE_EXTRA" = "y" ]; then

    ask_extra_args \
        "pbmm2 (Align)" \
        "pbmm2 align <reference> <input> <sample>.aligned.bam -j <cpus> --sort"
    ALIGN_ARGS="$EXTRA_ARGS_RETURN"

    ask_extra_args \
        "DeepVariant (CallVariants)" \
        "run_deepvariant --model_type PACBIO --ref <ref> --reads <bam> --output_vcf <sample>.variants.vcf.gz --intermediate_results_dir tmp --num_shards <cpus> --regions <regions>"
    DEEPVARIANT_ARGS="$EXTRA_ARGS_RETURN"

    if [ -n "$VEP_CACHE" ]; then
        ask_extra_args \
            "VEP (AnnotateWithVEP - cache mode)" \
            "vep --input_file <vcf> --output_file <out> --vcf --species ${VEP_SPECIES} --assembly ${VEP_ASSEMBLY} --cache --dir_cache <cache> --fasta <ref> --offline --fork <cpus> --sift b --symbol --gene_phenotype --regulatory --canonical --biotype --force_overwrite"
    else
        ask_extra_args \
            "VEP (AnnotateWithVEP - GFF mode)" \
            "vep --input_file <vcf> --output_file <out> --vcf --species ${VEP_SPECIES} --assembly ${VEP_ASSEMBLY} --gff <prepared.gff.gz> --fasta <ref> --fork <cpus> --symbol --canonical --biotype --force_overwrite"
    fi
    VEP_ARGS="$EXTRA_ARGS_RETURN"
fi

# =============================================================================
# Run options
# =============================================================================

read -p "Profile methylation with pb-CpG-tools (CpG + hap1/hap2)? (y/n) [n]: " METHYL_ANS
METHYL_ANS=${METHYL_ANS:-n}
read -p "Call structural variants with sawfish (DEL/INS/DUP/INV/BND)? (y/n) [n]: " SV_ANS
SV_ANS=${SV_ANS:-n}
read -p "Generate a global MultiQC report? (y/n) [y]: " MULTIQC_ANS
MULTIQC_ANS=${MULTIQC_ANS:-y}
read -p "Generate Nextflow execution reports (time/CPU/RAM per process)? (y/n) [y]: " EXEC_REPORTS
EXEC_REPORTS=${EXEC_REPORTS:-y}
read -p "Resume? (y/n) [n]: " RESUME
read -p "Clean Nextflow work cache after completion? (y/n) [n]: " CLEAN_CACHE

# =============================================================================
# Execution environment + resource ceiling.
# These rarely change between runs, so they are PERSISTED in run_profile.conf
# (next to run.sh) and reused on the next launch. Delete that file, or answer
# "reconfigure", to change them. Per-run inputs (sample, reference, intervals,
# run name, VEP source) are NOT stored here - they stay interactive.
# =============================================================================

PROFILE_FILE="${SCRIPT_DIR}/run_profile.conf"

prompt_execution_profile() {
    echo ""
    echo "=== Execution environment ==="
    echo "  Container engine and compute site are independent, composable profiles."
    echo "  Engine: singularity (clusters) | docker (workstation) | conda"
    echo "  Site:   hpc2n (SLURM) | local (run on this machine)"
    echo ""
    read -p "Container engine (singularity/docker/conda) [singularity]: " ENGINE
    ENGINE=${ENGINE:-singularity}
    read -p "Compute site (hpc2n/local): " SITE
    [ -z "$SITE" ] && echo "ERROR: Compute site required (hpc2n or local)" && exit 1

    SLURM_ACCOUNT=""
    if [ "$SITE" = "hpc2n" ]; then
        echo "SLURM account for '-A' (e.g. proj2026-1-23). Leave empty to use your"
        echo "cluster's default account."
        read -p "SLURM account: " SLURM_ACCOUNT
    fi

    SINGULARITY_CACHEDIR=""
    if [ "$ENGINE" = "singularity" ]; then
        echo "Shared Singularity image cache (downloaded once, reused across runs)."
        echo "Leave empty to use a repo-local cache (\${projectDir}/.singularity_cache)."
        read -p "Singularity cache dir [${NXF_SINGULARITY_CACHEDIR:-none}]: " SINGULARITY_CACHEDIR
        SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-${NXF_SINGULARITY_CACHEDIR:-}}
    fi

    echo ""
    echo "=== Resource allocation (heavy jobs: CallVariants / CallSV) ==="
    echo "  Only the two heavy jobs scale; everything else is fixed. The tier sets a"
    echo "  FIXED cpu/RAM request per job (time is fixed too). Pick a tier your cluster"
    echo "  can satisfy: a SLURM job that asks for more cpu/RAM than any free node has"
    echo "  will sit in the queue until a big enough node frees up."
    echo ""
    echo "  Make sure your cluster has at least the cpu/RAM shown for the chosen tier:"
    echo "    1) Normal     CallVariants 36 cpu / 120 GB   CallSV 36 cpu / 160 GB"
    echo "    2) High       CallVariants 48 cpu / 160 GB   CallSV 48 cpu / 200 GB"
    echo "    3) Very High  CallVariants 60 cpu / 200 GB   CallSV 60 cpu / 240 GB"
    echo "    4) Custom     you set the cpu/RAM (min 36 cpu / 160 GB), applied to both"
    echo ""
    echo "  NOTE: High and Very High request bigger nodes -> the jobs may WAIT LONGER in"
    echo "        the cluster queue. On a busy cluster, Normal (smaller nodes) often"
    echo "        finishes sooner overall. Bump the tier only if big nodes are readily"
    echo "        available on your cluster."
    read -p "Resource tier (1=Normal / 2=High / 3=Very High / 4=Custom) [1]: " RES_CHOICE
    RES_CHOICE=${RES_CHOICE:-1}

    CUSTOM_CPUS=""
    CUSTOM_MEM=""
    case "$RES_CHOICE" in
        1) RES_TIER="normal" ;;
        2) RES_TIER="high" ;;
        3) RES_TIER="veryhigh" ;;
        4)
            RES_TIER="custom"
            read -p "Custom CPUs per heavy job (min 36): " CUSTOM_CPUS
            if ! [[ "$CUSTOM_CPUS" =~ ^[0-9]+$ ]] || [ "$CUSTOM_CPUS" -lt 36 ]; then
                echo "ERROR: CPUs must be an integer >= 36"
                exit 1
            fi
            read -p "Custom memory per heavy job in GB (min 160): " CUSTOM_MEM
            if ! [[ "$CUSTOM_MEM" =~ ^[0-9]+$ ]] || [ "$CUSTOM_MEM" -lt 160 ]; then
                echo "ERROR: memory must be an integer (GB) >= 160"
                exit 1
            fi
            ;;
        *)
            echo "ERROR: Invalid choice (expected 1-4)"
            exit 1
            ;;
    esac
}

# Map the chosen tier to the per-process cpu/memory passed to nextflow.config.
# CallVariants and CallSV keep distinct baselines (120 vs 160 GB at Normal); High
# and Very High add the same delta to both. Custom sets both heavy jobs equal.
# Sets CV_CPUS / CV_MEM / SV_CPUS / SV_MEM (memory in whole GB).
compute_resources() {
    case "$RES_TIER" in
        normal)   CV_CPUS=36; CV_MEM=120; SV_CPUS=36; SV_MEM=160 ;;
        high)     CV_CPUS=48; CV_MEM=160; SV_CPUS=48; SV_MEM=200 ;;
        veryhigh) CV_CPUS=60; CV_MEM=200; SV_CPUS=60; SV_MEM=240 ;;
        custom)   CV_CPUS="$CUSTOM_CPUS"; CV_MEM="$CUSTOM_MEM"
                  SV_CPUS="$CUSTOM_CPUS"; SV_MEM="$CUSTOM_MEM" ;;
        *) echo "ERROR: unknown resource tier '$RES_TIER'" && exit 1 ;;
    esac
}

save_execution_profile() {
    cat > "$PROFILE_FILE" <<EOF
# Saved execution profile for run.sh - delete this file to reset.
# Edit by hand or answer "reconfigure" at the next run.sh launch.
ENGINE="$ENGINE"
SITE="$SITE"
SLURM_ACCOUNT="$SLURM_ACCOUNT"
SINGULARITY_CACHEDIR="$SINGULARITY_CACHEDIR"
RES_TIER="$RES_TIER"
CUSTOM_CPUS="$CUSTOM_CPUS"
CUSTOM_MEM="$CUSTOM_MEM"
EOF
    echo "Saved execution profile to: $PROFILE_FILE"
}

ENGINE="singularity"; SITE=""; SLURM_ACCOUNT=""; SINGULARITY_CACHEDIR=""
RES_TIER="normal"; CUSTOM_CPUS=""; CUSTOM_MEM=""

if [ -f "$PROFILE_FILE" ]; then
    # shellcheck disable=SC1090
    source "$PROFILE_FILE"
    echo ""
    echo "=== Execution profile (loaded from run_profile.conf) ==="
    echo "  Engine:        ${ENGINE}"
    echo "  Site:          ${SITE}"
    [ -n "${SLURM_ACCOUNT}" ]        && echo "  SLURM account: ${SLURM_ACCOUNT}"
    [ -n "${SINGULARITY_CACHEDIR}" ] && echo "  Sing. cache:   ${SINGULARITY_CACHEDIR}"
    if [ "${RES_TIER}" = "custom" ]; then
        echo "  Resource tier: custom (${CUSTOM_CPUS} cpu / ${CUSTOM_MEM} GB per heavy job)"
    else
        echo "  Resource tier: ${RES_TIER}"
    fi
    echo ""
    read -p "Reconfigure execution environment? (y/n) [n]: " RECONF
    if [ "${RECONF:-n}" = "y" ]; then
        prompt_execution_profile
        save_execution_profile
    fi
else
    prompt_execution_profile
    save_execution_profile
fi

# Guard against a partial/hand-edited run_profile.conf missing the engine or site.
[ -z "$ENGINE" ] && echo "ERROR: container engine missing (check run_profile.conf)" && exit 1
[ -z "$SITE" ]   && echo "ERROR: compute site missing (check run_profile.conf)" && exit 1
PROFILE="${ENGINE},${SITE}"

# Resolve the chosen tier into concrete cpu/memory for CallVariants and CallSV.
compute_resources

# Export the cache so nextflow.config's System.getenv('NXF_SINGULARITY_CACHEDIR') sees it.
if [ -n "${SINGULARITY_CACHEDIR}" ]; then
    export NXF_SINGULARITY_CACHEDIR="$SINGULARITY_CACHEDIR"
fi

# =============================================================================
# Resolve absolute paths
# =============================================================================

REFERENCE="$(readlink -f "$(eval echo "$REFERENCE")")"
[ -n "$VEP_CACHE" ] && VEP_CACHE="$(readlink -f "$(eval echo "$VEP_CACHE")")"
[ -n "$VEP_GFF" ]   && VEP_GFF="$(readlink -f "$(eval echo "$VEP_GFF")")"

mkdir -p "$RUN_DIR"

if [ "$MODE" = "single" ]; then
    INPUT_BAM="$(readlink -f "$(eval echo "$INPUT_BAM")")"
    # Generate a one-row-per-interval samplesheet for this single sample
    SAMPLESHEET="${RUN_DIR}/samplesheet.csv"
    echo "sample,bam,chromosome,interval" > "$SAMPLESHEET"
    for i in $(seq 0 $((NUM_INTERVALS - 1))); do
        echo "${SAMPLE_NAME},${INPUT_BAM},${CHROMS[$i]},${INTERVALS[$i]}" >> "$SAMPLESHEET"
    done
else
    SAMPLESHEET="$(readlink -f "$(eval echo "$SAMPLESHEET")")"
fi

# =============================================================================
# Summary
# =============================================================================

echo ""
echo "-- Summary --"
echo "Mode:         $MODE"
echo "Samplesheet:  $SAMPLESHEET"
if [ "$MODE" = "single" ]; then
    echo "Sample:       $SAMPLE_NAME"
    echo "Input BAM:    $INPUT_BAM"
fi
echo "Reference:    $REFERENCE"
if [ -n "$VEP_CACHE" ]; then
    echo "VEP source:   cache -> $VEP_CACHE"
else
    echo "VEP source:   GFF   -> $VEP_GFF"
fi
echo "Species:      $VEP_SPECIES"
echo "Assembly:     $VEP_ASSEMBLY"
echo "Run name:     $RUN_NAME"
echo "Output:       $RUN_DIR"
echo "Profile:      $PROFILE"
[ -n "$SLURM_ACCOUNT" ] && echo "SLURM acct:   $SLURM_ACCOUNT"
[ -n "$SINGULARITY_CACHEDIR" ] && echo "Sing. cache:  $SINGULARITY_CACHEDIR"
echo "Resources:    tier=${RES_TIER} | CallVariants ${CV_CPUS} cpu/${CV_MEM} GB | CallSV ${SV_CPUS} cpu/${SV_MEM} GB"
echo "Phasing:      HiPhase (automatic - phased VCF feeds VEP, haplotagged BAM produced)"
echo "Methylation:  $([ "$METHYL_ANS" = "y" ] && echo "pb-CpG-tools (combined + hap1/hap2, sliced + in report)" || echo "disabled")"
echo "SV calling:   $([ "$SV_ANS" = "y" ] && echo "sawfish (overlap-sliced + in report)" || echo "disabled")"
echo "MultiQC:      $([ "$MULTIQC_ANS" = "y" ] && echo "global report (${RUN_DIR}/multiqc/)" || echo "disabled")"
echo "Exec reports: $([ "$EXEC_REPORTS" = "y" ] && echo "${RUN_DIR}/pipeline_info/ (report+timeline+trace+dag)" || echo "disabled")"
echo ""
echo "Note: Variant calling regions are auto-detected from the reference genome."
echo "      Supported formats: Ensembl (1, 2, X...) and NCBI (NC_* accessions)."

if [ "$MODE" = "single" ]; then
    echo ""
    echo "Intervals to analyze:"
    for i in $(seq 0 $((NUM_INTERVALS - 1))); do
        echo "  $((i + 1)). chr${CHROMS[$i]}:${INTERVALS[$i]}"
    done
fi

if [ -n "${ALIGN_ARGS}${DEEPVARIANT_ARGS}${VEP_ARGS}" ]; then
    echo ""
    echo "Extra args:"
    [ -n "$ALIGN_ARGS" ]       && echo "  pbmm2:       $ALIGN_ARGS"
    [ -n "$DEEPVARIANT_ARGS" ] && echo "  DeepVariant: $DEEPVARIANT_ARGS"
    [ -n "$VEP_ARGS" ]         && echo "  VEP:         $VEP_ARGS"
fi

echo ""
read -p "Run pipeline? (y/n): " CONFIRM
[ "$CONFIRM" != "y" ] && echo "Aborted." && exit 0

# =============================================================================
# Build Nextflow arguments
# =============================================================================

ARGS=(
    -profile       "$PROFILE"
    --input        "$SAMPLESHEET"
    --reference    "$REFERENCE"
    --outdir       "$OUTPUT_DIR"
    --run_name     "$RUN_NAME"
    --vep_species  "$VEP_SPECIES"
    --vep_assembly "$VEP_ASSEMBLY"
    --callvariants_cpus   "$CV_CPUS"
    --callvariants_memory "${CV_MEM}.GB"
    --callsv_cpus         "$SV_CPUS"
    --callsv_memory       "${SV_MEM}.GB"
)

[ -n "$SLURM_ACCOUNT" ] && ARGS+=(--slurm_account "$SLURM_ACCOUNT")

[ -n "$VEP_CACHE" ] && ARGS+=(--vep_cache "$VEP_CACHE")
[ -n "$VEP_GFF" ]   && ARGS+=(--vep_gff   "$VEP_GFF")

[ "$METHYL_ANS"  = "y" ] && ARGS+=(--run_methylation true)
[ "$SV_ANS"      = "y" ] && ARGS+=(--run_sv true)
[ "$MULTIQC_ANS" != "y" ] && ARGS+=(--skip_multiqc true)

# Nextflow execution reports (time/CPU/RAM per process). Timestamped so repeated
# runs don't clash (Nextflow refuses to overwrite an existing report file).
if [ "$EXEC_REPORTS" = "y" ]; then
    REPORT_DIR="${RUN_DIR}/pipeline_info"
    mkdir -p "$REPORT_DIR"
    STAMP="$(date +%Y%m%d_%H%M%S)"
    ARGS+=(
        -with-report   "${REPORT_DIR}/report_${STAMP}.html"
        -with-timeline "${REPORT_DIR}/timeline_${STAMP}.html"
        -with-trace    "${REPORT_DIR}/trace_${STAMP}.txt"
        -with-dag      "${REPORT_DIR}/dag_${STAMP}.html"
    )
fi

[ "$RESUME" = "y" ] && ARGS+=(-resume)

# Pass extra args (quoted as a single string each)
[ -n "$ALIGN_ARGS" ]       && ARGS+=(--align_args       "$ALIGN_ARGS")
[ -n "$DEEPVARIANT_ARGS" ] && ARGS+=(--deepvariant_args "$DEEPVARIANT_ARGS")
[ -n "$VEP_ARGS" ]         && ARGS+=(--vep_args         "$VEP_ARGS")

# =============================================================================
# Run the pipeline
# =============================================================================

nextflow run "${SCRIPT_DIR}/find_specific_mutation.nf" "${ARGS[@]}"

# Point the user to the run-level reports (aggregated summary + MultiQC at the run
# root, and the per-interval HTML reports inside each sample folder).
echo ""
echo "Reports:"
echo "  - Aggregated report: ${RUN_DIR}/general_summary.html"
echo "  - MultiQC report:    ${RUN_DIR}/multiqc_report.html"
echo "  - Per-interval:      ${RUN_DIR}/<sample>/"

# =============================================================================
# Post-run cleanup (Nextflow work cache only - does NOT touch VEP cache or GFF)
# =============================================================================

if [ "$CLEAN_CACHE" = "y" ]; then
    echo ""
    echo "Cleaning Nextflow work cache for the last run..."
    nextflow clean -f last
    echo "Done."
fi
