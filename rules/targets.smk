# -------------------------------------------------------------------------
# Final outputs
# -------------------------------------------------------------------------
rule all:
    default_target: True
    input:
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

        # Mapping and depth
        *(
            expand(alignment_path("{sample}", "{sample}.vs_assembly.sorted.bam"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED else []
        ),
        *(
            expand(alignment_path("{sample}", "{sample}.vs_assembly.sorted.bam.bai"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED else []
        ),
        *(
            expand(alignment_path("{sample}", "{sample}.depth.txt"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED else []
        ),

        # Binning
        *(
            expand(binning_path("MetaBAT2", "{sample}", ".done"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "metabat2" in BINNING_TOOLS else []
        ),
        *(
            expand(binning_path("SemiBin2", "{sample}", ".done"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "semibin2" in BINNING_TOOLS else []
        ),
        *(
            expand(binning_path("COMEBin", "{sample}", ".done"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "comebin" in BINNING_TOOLS else []
        ),

        # Collected bins
        *(
            expand(binning_path("AllBins", "{sample}", ".done"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and len(BINNING_TOOLS) > 0 else []
        ),

        # Normalized contig-to-bin tables
        *(
            expand(binning_path("Normalized", "{sample}", "metabat2.contig2bin.tsv"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "metabat2" in BINNING_TOOLS else []
        ),
        *(
            expand(binning_path("Normalized", "{sample}", "semibin2.contig2bin.tsv"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "semibin2" in BINNING_TOOLS else []
        ),
        *(
            expand(binning_path("Normalized", "{sample}", "comebin.contig2bin.tsv"), sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "comebin" in BINNING_TOOLS else []
        ),

        # Prokaryotic refinement
        *(
            expand(refinement_path("Prokaryotic", "DASTool", "{sample}", "{sample}_DASTool_summary.tsv"), sample=SAMPLES)
            if DASTOOL_ENABLED else []
        ),
        *(
            expand(refinement_path("Prokaryotic", "DASTool", "{sample}", "{sample}_DASTool_contig2bin.tsv"), sample=SAMPLES)
            if DASTOOL_ENABLED else []
        ),
        *(
            expand(refinement_path("Prokaryotic", "DASTool", "{sample}", "{sample}_DASTool_bins"), sample=SAMPLES)
            if DASTOOL_ENABLED else []
        ),
        *(
            expand(refinement_path("Prokaryotic", "DASTool", "{sample}", ".done"), sample=SAMPLES)
            if DASTOOL_ENABLED else []
        ),

        # Prokaryotic quality estimation
        *(
            expand(refinement_path("Prokaryotic", "CheckM2", "{sample}", "quality_report.tsv"), sample=SAMPLES)
            if CHECKM2_ENABLED else []
        ),
        *(
            expand(refinement_path("Prokaryotic", "CheckM2", "{sample}", ".done"), sample=SAMPLES)
            if CHECKM2_ENABLED else []
        ),

        # Prokaryotic classification
        *(
            expand(refinement_path("Prokaryotic", "GTDBTk", "{sample}", "gtdbtk.bac120.summary.tsv"), sample=SAMPLES)
            if GTDBTK_ENABLED else []
        ),
        *(
            expand(refinement_path("Prokaryotic", "GTDBTk", "{sample}", "gtdbtk.ar53.summary.tsv"), sample=SAMPLES)
            if GTDBTK_ENABLED else []
        ),
        *(
            expand(refinement_path("Prokaryotic", "GTDBTk", "{sample}", ".done"), sample=SAMPLES)
            if GTDBTK_ENABLED else []
        ),

        # Eukaryotic refinement
        *(
            expand(refinement_path("Eukaryotic", "ACR", "{sample}", ".done"), sample=SAMPLES)
            if ACR_ENABLED else []
        ),
        *(
            expand(refinement_path("Eukaryotic", "EukBins", "{sample}", "kept_bins.tsv"), sample=SAMPLES)
            if ACR_ENABLED else []
        ),
        *(
            expand(refinement_path("Eukaryotic", "EukBins", "{sample}", ".done"), sample=SAMPLES)
            if ACR_ENABLED else []
        ),

        # Eukaryotic dereplication
        *(
            expand(refinement_path("Eukaryotic", "dRep", "{sample}", "data_tables", "Cdb.csv"), sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED else []
        ),
        *(
            expand(refinement_path("Eukaryotic", "dRep", "{sample}", ".done"), sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED else []
        ),

        # Eukaryotic quality estimation
        *(
            expand(refinement_path("Eukaryotic", "EukCC", "{sample}", "eukcc.csv"), sample=SAMPLES)
            if ACR_ENABLED and EUKCC_ENABLED else []
        ),
        *(
            expand(refinement_path("Eukaryotic", "EukCC", "{sample}", ".done"), sample=SAMPLES)
            if ACR_ENABLED and EUKCC_ENABLED else []
        ),

        # Final eukaryotic bin selection
        *(
            expand(refinement_path("Eukaryotic", "FinalBins", "{sample}", "selected_bins.tsv"), sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED else []
        ),
        *(
            expand(refinement_path("Eukaryotic", "FinalBins", "{sample}", ".done"), sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED else []
        ),

        # Eukaryotic bin classification
        *(
            expand(refinement_path("Eukaryotic", "BAT", "{sample}", "bin2classification.txt"), sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED and BAT_ENABLED else []
        ),
        *(
            expand(refinement_path("Eukaryotic", "BAT", "{sample}", ".done"), sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED and BAT_ENABLED else []
        ),

        # Genome-level final reports
        *(
            expand(report_path("{sample}", "{sample}.genome_inventory.tsv"), sample=SAMPLES)
            if RUN_ASSEMBLY else []
        ),
        *(
            expand(report_path("{sample}", "{sample}.genome_summary.tsv"), sample=SAMPLES)
            if RUN_ASSEMBLY else []
        ),
        *(
            expand(report_path("{sample}", "{sample}.bin_trace.tsv"), sample=SAMPLES)
            if RUN_ASSEMBLY else []
        ),

        # Combined pipeline summary
        *(
            expand(report_path("{sample}", "{sample}.pipeline_summary.tsv"), sample=SAMPLES)
            if RUN_PROFILE and RUN_ASSEMBLY else []
        ),

        # HTML report (generated for all active pipeline modes)
        *(
            expand(report_path("{sample}", "{sample}.taxonomy_top10.json"), sample=SAMPLES)
            if TAXONOMY_ENABLED else []
        ),
        os.path.join("Reports", "report.html")
