# -------------------------------------------------------------------------
# Final outputs
# -------------------------------------------------------------------------
# Wrapped in a function (rather than a static list) because some entries
# depend on the assembly_binning_gate checkpoint via binning_samples(),
# which can only be evaluated lazily once assembly has completed.
def all_targets(wildcards):
    return [
        # -----------------------------------------------------------------
        # Versions manifest (written to Reports/ by run wrapper)
        # -----------------------------------------------------------------
        os.path.join("Reports", "versions.json"),

        # -----------------------------------------------------------------
        # Raw QC
        # -----------------------------------------------------------------
        *(
            expand(qc_path("Raw", "{sample}", "{sample}_NanoPlot-report.html"), sample=SAMPLES)
            if QC_TOOL == "nanoplot" else []
        ),
        *(
            expand(qc_path("Raw", "{sample}", "{sample}_nanoQC.html"), sample=SAMPLES)
            if QC_TOOL == "nanoqc" else []
        ),

        # -----------------------------------------------------------------
        # Preprocessing final
        # -----------------------------------------------------------------
        *(
            expand(
                preprocessing_path("{sample}", "{sample}." + FILTERING_METHOD + ".fastq.gz"),
                sample=SAMPLES
            )
            if PREPROCESSING_ENABLED else []
        ),
        # -----------------------------------------------------------------
        # Host removal stats
        # -----------------------------------------------------------------
        *(
            expand(preprocessing_path("{sample}", "{sample}.host_removal_stats.json"), sample=SAMPLES)
            if HOST_REMOVAL_ENABLED else []
        ),

        # -----------------------------------------------------------------
        # Filtered QC
        # -----------------------------------------------------------------
        *(
            expand(qc_path("Filtered", "{sample}", "{sample}_NanoPlot-report.html"), sample=SAMPLES)
            if PREPROCESSING_ENABLED and QC_TOOL == "nanoplot" else []
        ),
        *(
            expand(qc_path("Filtered", "{sample}", "{sample}_nanoQC.html"), sample=SAMPLES)
            if PREPROCESSING_ENABLED and QC_TOOL == "nanoqc" else []
        ),

        # -----------------------------------------------------------------
        # Profiling branch
        # -----------------------------------------------------------------
        *(
            expand(taxonomy_path("Kraken2", "{sample}", "{sample}.kraken2.output.txt"), sample=SAMPLES)
            if TAXONOMY_ENABLED and "kraken2" in TAX_TOOLS else []
        ),
        *(
            expand(taxonomy_path("Kraken2", "{sample}", "{sample}.kraken2.report.txt"), sample=SAMPLES)
            if TAXONOMY_ENABLED and "kraken2" in TAX_TOOLS else []
        ),
        *(
            expand(taxonomy_path("Centrifuger", "{sample}", "{sample}.centrifuger.classification.tsv"), sample=SAMPLES)
            if TAXONOMY_ENABLED and "centrifuger" in TAX_TOOLS else []
        ),
        *(
            expand(taxonomy_path("Centrifuger", "{sample}", "{sample}.centrifuger.report.tsv"), sample=SAMPLES)
            if TAXONOMY_ENABLED and "centrifuger" in TAX_TOOLS else []
        ),
        *(
            expand(stage_path("Profile", "{sample}", "{sample}.done"), sample=SAMPLES)
            if RUN_PROFILE and RUN_ASSEMBLY and SERIAL_PROFILE_THEN_ASSEMBLY else []
        ),

        # -----------------------------------------------------------------
        # Assembly / binning / genome branch
        # -----------------------------------------------------------------
        *(
            expand(assembly_path("MetaFlye", "{sample}", "assembly.fasta"), sample=SAMPLES)
            if ASSEMBLY_ENABLED else []
        ),
        *(
            expand(assembly_path("MetaFlye", "{sample}", "assembly_info.txt"), sample=SAMPLES)
            if ASSEMBLY_ENABLED else []
        ),
        *(
            expand(assembly_path("MetaFlye", "{sample}", "assembly_graph.gfa"), sample=SAMPLES)
            if ASSEMBLY_ENABLED else []
        ),

        # Assembly QC
        *(
            expand(assembly_path("MetaQUAST", "{sample}", "report.html"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and METAQUAST_ENABLED else []
        ),
        *(
            expand(assembly_path("MetaQUAST", "{sample}", "report.tsv"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and METAQUAST_ENABLED else []
        ),

        # Mapping and depth (only for samples whose assembly passes the
        # post-assembly binning gate; mapping/depth are only used for binning)
        *(
            [alignment_path(s, f"{s}.vs_assembly.sorted.bam") for s in binning_samples()]
            if ASSEMBLY_ENABLED and MAPPING_ENABLED else []
        ),
        *(
            [alignment_path(s, f"{s}.vs_assembly.sorted.bam.bai") for s in binning_samples()]
            if ASSEMBLY_ENABLED and MAPPING_ENABLED else []
        ),
        *(
            [alignment_path(s, f"{s}.depth.txt") for s in binning_samples()]
            if ASSEMBLY_ENABLED and MAPPING_ENABLED else []
        ),

        # Binning (skipped for samples whose assembly does not have enough
        # contigs of sufficient length, per the post-assembly binning gate)
        *(
            [binning_path("MetaBAT2", s, ".done") for s in binning_samples()]
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "metabat2" in BINNING_TOOLS else []
        ),
        *(
            [binning_path("SemiBin2", s, ".done") for s in binning_samples()]
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "semibin2" in BINNING_TOOLS else []
        ),
        *(
            [binning_path("COMEBin", s, ".done") for s in binning_samples()]
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "comebin" in BINNING_TOOLS else []
        ),

        # Collected bins
        *(
            [binning_path("AllBins", s, ".done") for s in binning_samples()]
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and len(BINNING_TOOLS) > 0 else []
        ),

        # Normalized contig-to-bin tables
        *(
            [binning_path("Normalized", s, "metabat2.contig2bin.tsv") for s in binning_samples()]
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "metabat2" in BINNING_TOOLS else []
        ),
        *(
            [binning_path("Normalized", s, "semibin2.contig2bin.tsv") for s in binning_samples()]
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "semibin2" in BINNING_TOOLS else []
        ),
        *(
            [binning_path("Normalized", s, "comebin.contig2bin.tsv") for s in binning_samples()]
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "comebin" in BINNING_TOOLS else []
        ),

        # Prokaryotic refinement (skipped for samples that did not pass the
        # post-assembly binning gate, since binning outputs are unavailable)
        *(
            [refinement_path("Prokaryotic", "DASTool", s, f"{s}_DASTool_summary.tsv") for s in binning_samples()]
            if DASTOOL_ENABLED else []
        ),
        *(
            [refinement_path("Prokaryotic", "DASTool", s, f"{s}_DASTool_contig2bin.tsv") for s in binning_samples()]
            if DASTOOL_ENABLED else []
        ),
        *(
            [refinement_path("Prokaryotic", "DASTool", s, f"{s}_DASTool_bins") for s in binning_samples()]
            if DASTOOL_ENABLED else []
        ),
        *(
            [refinement_path("Prokaryotic", "DASTool", s, ".done") for s in binning_samples()]
            if DASTOOL_ENABLED else []
        ),

        # Prokaryotic quality estimation
        *(
            [refinement_path("Prokaryotic", "CheckM2", s, "quality_report.tsv") for s in binning_samples()]
            if CHECKM2_ENABLED else []
        ),
        *(
            [refinement_path("Prokaryotic", "CheckM2", s, ".done") for s in binning_samples()]
            if CHECKM2_ENABLED else []
        ),

        # Prokaryotic classification
        *(
            [refinement_path("Prokaryotic", "GTDBTk", s, "gtdbtk.bac120.summary.tsv") for s in binning_samples()]
            if GTDBTK_ENABLED else []
        ),
        *(
            [refinement_path("Prokaryotic", "GTDBTk", s, "gtdbtk.ar53.summary.tsv") for s in binning_samples()]
            if GTDBTK_ENABLED else []
        ),
        *(
            [refinement_path("Prokaryotic", "GTDBTk", s, ".done") for s in binning_samples()]
            if GTDBTK_ENABLED else []
        ),

        # Eukaryotic refinement
        *(
            [refinement_path("Eukaryotic", "ACR", s, ".done") for s in binning_samples()]
            if ACR_ENABLED else []
        ),
        *(
            [refinement_path("Eukaryotic", "EukBins", s, "kept_bins.tsv") for s in binning_samples()]
            if ACR_ENABLED else []
        ),
        *(
            [refinement_path("Eukaryotic", "EukBins", s, ".done") for s in binning_samples()]
            if ACR_ENABLED else []
        ),

        # Eukaryotic dereplication
        *(
            [refinement_path("Eukaryotic", "dRep", s, "data_tables", "Cdb.csv") for s in binning_samples()]
            if ACR_ENABLED and DREP_ENABLED else []
        ),
        *(
            [refinement_path("Eukaryotic", "dRep", s, ".done") for s in binning_samples()]
            if ACR_ENABLED and DREP_ENABLED else []
        ),

        # Eukaryotic quality estimation
        *(
            [refinement_path("Eukaryotic", "EukCC", s, "eukcc.csv") for s in binning_samples()]
            if ACR_ENABLED and EUKCC_ENABLED else []
        ),
        *(
            [refinement_path("Eukaryotic", "EukCC", s, ".done") for s in binning_samples()]
            if ACR_ENABLED and EUKCC_ENABLED else []
        ),

        # Final eukaryotic bin selection
        *(
            [refinement_path("Eukaryotic", "FinalBins", s, "selected_bins.tsv") for s in binning_samples()]
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED else []
        ),
        *(
            [refinement_path("Eukaryotic", "FinalBins", s, ".done") for s in binning_samples()]
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED else []
        ),

        # Eukaryotic bin classification
        *(
            [refinement_path("Eukaryotic", "BAT", s, "bin2classification.txt") for s in binning_samples()]
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED and BAT_ENABLED else []
        ),
        *(
            [refinement_path("Eukaryotic", "BAT", s, ".done") for s in binning_samples()]
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED and BAT_ENABLED else []
        ),

        # Genome-level final reports (only for samples that passed the
        # binning gate; a skipped sample has no bins to report on)
        *(
            [report_path(s, f"{s}.genome_inventory.tsv") for s in binning_samples()]
            if RUN_ASSEMBLY else []
        ),
        *(
            [report_path(s, f"{s}.genome_summary.tsv") for s in binning_samples()]
            if RUN_ASSEMBLY else []
        ),
        *(
            [report_path(s, f"{s}.bin_trace.tsv") for s in binning_samples()]
            if RUN_ASSEMBLY else []
        ),

        # Combined pipeline summary
        *(
            [report_path(s, f"{s}.pipeline_summary.tsv") for s in binning_samples()]
            if RUN_PROFILE and RUN_ASSEMBLY else []
        ),

        # HTML report (generated for all active pipeline modes)
        *(
            expand(report_path("{sample}", "{sample}.taxonomy_top10.json"), sample=SAMPLES)
            if TAXONOMY_ENABLED else []
        ),
        os.path.join("Reports", "report.html")
    ]


rule all:
    default_target: True
    input:
        all_targets
