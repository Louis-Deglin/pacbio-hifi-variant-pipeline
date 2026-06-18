include { DetectAndStandardize } from './modules/detect_and_standardize.nf'
include { Align } from './modules/alignment.nf'
include { CallVariants } from './modules/variant_calling.nf'
include { HiPhase } from './modules/phasing.nf'
include { ExtractIntervals } from './modules/extract_intervals.nf'
include { MethylationProfiling } from './modules/methylation.nf'
include { SliceMethylation } from './modules/slice_methylation.nf'
include { CallSV } from './modules/call_sv.nf'
include { SliceSV } from './modules/slice_sv.nf'
include { AnnotateSVWithVEP } from './modules/annotate_sv_vep.nf'
include { IntervalReport } from './modules/interval_report.nf'
include { RunReport } from './modules/run_report.nf'
include { AnnotateWithVEP } from './modules/annotate_vep.nf'
include { SamtoolsStats } from './modules/samtools_stats.nf'
include { BcftoolsStats } from './modules/bcftools_stats.nf'
include { BcftoolsStats as BcftoolsStatsSV } from './modules/bcftools_stats.nf'
include { MultiQC } from './modules/multiqc.nf'

workflow {
    main:
    // Parse + validate the CSV samplesheet.
    // Columns: sample,bam,chromosome,interval - one row per (sample, interval).
    // Several intervals for the same sample = repeat the row with the same
    // sample + bam and a different chromosome/interval.
    def samples_ch = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample     = row.sample?.trim()
            def bam        = row.bam?.trim()
            def chromosome = row.chromosome?.trim()
            def interval   = row.interval?.trim()
            if (!sample || !bam || !chromosome || !interval) {
                error "Samplesheet row with an empty field (need sample,bam,chromosome,interval): ${row}"
            }
            // file() does not expand '~' - BAM paths must be absolute.
            tuple(sample, file(bam, checkIfExists: true), chromosome, interval)
        }

    reference_ch = channel.fromPath(params.reference, checkIfExists: true)

    // Detect format (Ensembl/NCBI), standardize FASTA, index, extract regions - runs ONCE.
    DetectAndStandardize(reference_ch)

    // Shared reference is reused by every sample -> value channels (.first()),
    // otherwise the single-item queue channel is consumed by the first sample.
    def std_value      = DetectAndStandardize.out.standardized.first()
    def regions_value  = DetectAndStandardize.out.regions.first()
    def ref_only_value = std_value.map { ref, fai -> ref }

    // One alignment + variant call per unique (sample, bam), even when the sample
    // appears on several rows (multiple intervals).
    def reads_ch = samples_ch
        .map { sample, bam, chromosome, interval -> tuple(sample, bam) }
        .unique()

    // Align reads to the standardized reference (pbmm2 only needs the .fa)
    Align(reads_ch, ref_only_value)

    // Variant calling with aligned BAM, standardized reference+index, detected regions
    CallVariants(Align.out.bam, std_value, regions_value)

    // Phasing : HiPhase needs the aligned BAM + the called VCF, joined by
    // sample. It produces a phased VCF (GT 0|1) used downstream AND a haplotagged BAM
    // (HP/PS tags) reused by the optional methylation rail.
    def phase_input = Align.out.bam.join(CallVariants.out.vcf)
    HiPhase(phase_input, std_value)

    // Intervals stay attached to their sample; join each sample's PHASED vcf back to
    // its interval(s). combine(by: 0) gives one (sample, vcf, tbi, chrom, interval) per interval.
    def intervals_ch = samples_ch
        .map { sample, bam, chromosome, interval -> tuple(sample, chromosome, interval) }

    def vcf_with_intervals = HiPhase.out.phased_vcf.combine(intervals_ch, by: 0)
    ExtractIntervals(vcf_with_intervals)

    // Choose VEP annotation source: cache directory or GFF file (reused by all samples)
    def annotation_ch
    if (params.vep_cache) {
        annotation_ch = channel.fromPath(params.vep_cache, type: 'dir', checkIfExists: true).first()
    } else {
        annotation_ch = channel.fromPath(params.vep_gff, checkIfExists: true).first()
    }

    // VEP needs both the .fa and its .fai staged in the work dir
    def vep_input = ExtractIntervals.out.interval_vcf
        .combine(annotation_ch)
        .combine(std_value)

    AnnotateWithVEP(vep_input)

    // Methylation profiling (opt-in via --run_methylation). Branches off the HiPhase
    // haplotagged BAM so pb-CpG-tools can emit hap1/hap2 tracks (allele-specific 5mC).
    // The genome-wide BEDs are sliced per interval; SliceMethylation also produces a
    // per-interval summary (stats + CpG table) that is embedded in the IntervalReport HTML.
    def methylation_ch        = channel.empty()   // genome-wide BEDs (combined/hap1/hap2)
    def methylation_bw_ch     = channel.empty()   // bigwig tracks
    def methylation_sliced_ch = channel.empty()   // per-interval sliced BEDs
    def meth_summary_ch       = channel.empty()   // per-interval summary text

    if (params.run_methylation) {
        MethylationProfiling(HiPhase.out.haplotagged_bam, std_value)

        // Slice every track (combined/hap1/hap2) to each interval of the sample.
        def meth_with_intervals = MethylationProfiling.out.beds.combine(intervals_ch, by: 0)
        SliceMethylation(meth_with_intervals)

        methylation_ch        = MethylationProfiling.out.beds
        methylation_bw_ch     = MethylationProfiling.out.bigwig
        methylation_sliced_ch = SliceMethylation.out.sliced
        meth_summary_ch       = SliceMethylation.out.summary
    }

    // Structural variant calling (opt-in via --run_sv). Branches off the ALIGNED BAM
    // (independent of phasing) so sawfish sees the raw alignments. The genome-wide SV VCF
    // is sliced per interval with overlap semantics; SliceSV also produces a per-interval
    // summary (counts by SVTYPE + per-SV table) embedded in the IntervalReport HTML.
    def sv_ch            = channel.empty()   // genome-wide SV VCF
    def sv_annotated_ch  = channel.empty()   // per-interval VEP-annotated SV VCF
    def sv_summary_ch    = channel.empty()   // per-interval summary text (built from annotated VCF)
    def sv_vep_summary_ch = channel.empty()  // VEP *_summary.html for MultiQC

    if (params.run_sv) {
        CallSV(Align.out.bam, std_value)

        // Slice the SV VCF to each interval of the sample (overlap, not containment).
        def sv_with_intervals = CallSV.out.sv_vcf.combine(intervals_ch, by: 0)
        SliceSV(sv_with_intervals)

        // Annotate the sliced SV VCF with VEP (same gene model as the SNV path). VEP reads
        // symbolic alleles/BND and assigns SV consequence terms; AnnotateSVWithVEP also builds
        // the per-interval summary (with gene/consequence) embedded in the variants_list report.
        def sv_vep_input = SliceSV.out.sliced
            .combine(annotation_ch)
            .combine(std_value)
        AnnotateSVWithVEP(sv_vep_input)

        sv_ch             = CallSV.out.sv_vcf
        sv_annotated_ch   = AnnotateSVWithVEP.out.annotated_sv
        sv_summary_ch     = AnnotateSVWithVEP.out.summary
        sv_vep_summary_ch = AnnotateSVWithVEP.out.vep_summary
    }

    // Per-(sample x interval) HTML report. Attach the methylation summary AND the
    // annotated SV VCF by (sample, chromosome, interval); join(remainder: true) keeps every
    // annotated SNV VCF even when methylation/SV are off/absent, falling back to placeholders.
    // A path() input needs a real file, so the placeholders are generated in the Nextflow
    // work dir (not committed to the repo, cleaned by `nextflow clean`). The render script
    // keys off params.run_methylation/run_sv, so the placeholder file content is irrelevant
    // when a rail is off. The HTML parser reads the VCF with gzip, so no .tbi is needed.
    def no_meth = file("${workflow.workDir}/NO_METHYLATION.txt")
    no_meth.text = "Methylation: not available for this region (profiling disabled via --run_methylation, or no MM/ML tags in the BAM).\n"
    def no_sv = file("${workflow.workDir}/NO_SV.txt")
    no_sv.text = "Structural variants: not available for this region (SV calling disabled via --run_sv).\n"

    def vcf_keyed  = AnnotateWithVEP.out.annotated_vcf
        .map { s, c, i, vcf, tbi -> tuple([s, c, i], vcf) }
    def meth_keyed = meth_summary_ch
        .map { s, c, i, sum -> tuple([s, c, i], sum) }
    def sv_keyed = sv_annotated_ch
        .map { s, c, i, vcf, tbi -> tuple([s, c, i], vcf) }

    def ir_input = vcf_keyed
        .join(meth_keyed, remainder: true)
        .join(sv_keyed,   remainder: true)
        .map { key, snv_vcf, meth, sv_vcf -> tuple(key[0], key[1], key[2], snv_vcf, meth ?: no_meth, sv_vcf ?: no_sv) }

    def render_script = file("${projectDir}/bin/render_interval.py")
    IntervalReport(ir_input, render_script)

    // Run-level HTML (built ONCE): an aggregated general_summary.html over every interval.
    // The fragments are collected; a manifest (sample,chrom,interval) drives the plain-text
    // table of contents, reconstructing fragment filenames deterministically.
    def report_manifest = IntervalReport.out.report
        .map { s, c, i, html -> "${s}\t${c}\t${i}" }
        .collectFile(name: 'report_manifest.tsv', newLine: true, sort: true)
    def run_report_script = file("${projectDir}/bin/build_run_report.py")
    RunReport(IntervalReport.out.fragment.collect(), report_manifest, run_report_script)

    // QC + global MultiQC report (one report aggregating all samples).
    // Disable with --skip_multiqc.
    def multiqc_ch = channel.empty()
    if (!params.skip_multiqc) {
        SamtoolsStats(Align.out.bam)
        BcftoolsStats(CallVariants.out.vcf.map { s, vcf -> tuple(s, vcf, '') })

        // bcftools stats on the SV calls too (counts by type/size). Reuses BcftoolsStats via
        // an alias with a 'sv' tag (distinct filename). sv_ch is empty when --run_sv is off,
        // so the aliased process simply doesn't run - no extra guard needed.
        BcftoolsStatsSV(sv_ch.map { s, vcf -> tuple(s, vcf, 'sv') })

        // Gather every QC artefact (alignment stats, SNV + SV variant stats, VEP summaries)
        // into a single channel, then collect() so MultiQC runs ONCE over all samples.
        def qc_files = SamtoolsStats.out.stats
            .mix(BcftoolsStats.out.stats)
            .mix(BcftoolsStatsSV.out.stats)
            .mix(AnnotateWithVEP.out.vep_summary)
            .mix(sv_vep_summary_ch)
            .collect()

        MultiQC(qc_files)
        // Publish only the standalone HTML (at the run root); drop the multiqc_data/ dir.
        multiqc_ch = MultiQC.out.report.map { report, data -> report }
    }

    publish:
    aligned_bam        = Align.out.bam
    raw_vcf            = CallVariants.out.vcf
    deepvariant_report = CallVariants.out.report
    phasing            = HiPhase.out.phased_vcf
    phasing_stats      = HiPhase.out.stats
    haplotagged_bam    = HiPhase.out.haplotagged_bam
    intervals_vcf      = ExtractIntervals.out.interval_vcf
    vep_annotated      = AnnotateWithVEP.out.annotated_vcf
    report             = IntervalReport.out.report
    run_summary        = RunReport.out.summary
    methylation        = methylation_ch
    methylation_bw     = methylation_bw_ch
    methylation_sliced = methylation_sliced_ch
    sv                 = sv_ch
    sv_annotated       = sv_annotated_ch
    multiqc            = multiqc_ch
}

output {
    // Aligned BAM + index (pbmm2). Heavy, but published on purpose for a complete results/.
    aligned_bam {
        path { sample, bam, bai -> "${sample}/alignment" }
        mode 'copy'
    }

    // Raw DeepVariant VCF (pre-phasing, pre-interval).
    raw_vcf {
        path { sample, vcf -> "${sample}/variant_calling" }
        mode 'copy'
    }

    // DeepVariant standalone QC HTML (not a MultiQC input - MultiQC has no DeepVariant module).
    deepvariant_report {
        path { sample, html -> "${sample}/variant_calling" }
        mode 'copy'
    }

    phasing {
        path { sample, vcf, tbi -> "${sample}/phasing" }
        mode 'copy'
    }

    phasing_stats {
        path { sample, stats, summary, blocks -> "${sample}/phasing" }
        mode 'copy'
    }

    // HiPhase haplotagged BAM + index (HP/PS tags). Heavy, published for a complete results/.
    haplotagged_bam {
        path { sample, bam, bai -> "${sample}/phasing" }
        mode 'copy'
    }

    // Per-interval VCF sliced out before VEP annotation.
    intervals_vcf {
        path { sample, chromosome, interval, vcf, tbi -> "${sample}/intervals" }
        mode 'copy'
    }

    vep_annotated {
        path { sample, chromosome, interval, vcf, tbi -> "${sample}/vep_annotated" }
        mode 'copy'
    }

    // Per-(sample x interval) HTML report, published directly under the sample folder
    // (alongside alignment/, phasing/, ... - no extra report/ subfolder).
    report {
        path { sample, chromosome, interval, html -> "${sample}" }
        mode 'copy'
    }

    // Run-level aggregated report, at the run root (same level as multiqc_report.html).
    run_summary {
        path { html -> "." }
        mode 'copy'
    }

    // Methylation: genome-wide BEDs + bigwig, plus per-interval sliced BEDs.
    methylation {
        path { sample, beds -> "${sample}/methylation" }
        mode 'copy'
    }

    methylation_bw {
        path { sample, bw -> "${sample}/methylation" }
        mode 'copy'
    }

    methylation_sliced {
        path { sample, chromosome, interval, beds -> "${sample}/methylation" }
        mode 'copy'
    }

    // Structural variants: genome-wide SV VCF + per-interval VEP-annotated VCFs.
    sv {
        path { sample, vcf -> "${sample}/sv" }
        mode 'copy'
    }

    sv_annotated {
        path { sample, chromosome, interval, vcf, tbi -> "${sample}/sv" }
        mode 'copy'
    }

    // Single global report shared by all samples: just the standalone HTML at the run
    // root (the multiqc_data/ directory is intentionally not published).
    multiqc {
        path { report -> "." }
        mode 'copy'
    }
}
