# -------------------------------------------------------------------------
# Final outputs
# -------------------------------------------------------------------------
rule all:
    default_target: True
    input:
        # -----------------------------------------------------------------
        # Raw QC
        # -----------------------------------------------------------------
        expand("outputs/{sample}/reports/QC/{sample}_RAW_NanoPlot-report.html", sample=SAMPLES),
        expand("outputs/{sample}/reports/QC/{sample}_RAW_nanoQC.html", sample=SAMPLES),

        # -----------------------------------------------------------------
        # Preprocessing final
        # -----------------------------------------------------------------
        # Only requested when preprocessing is enabled in config.yaml:
        #
        # preprocessing:
        #   enabled: true
        *(
            expand("outputs/{sample}/reports/preprocessing/{sample}-preprocessed.fastq.gz", sample=SAMPLES)
            if PREPROCESSING_ENABLED else []
        ),

        # -----------------------------------------------------------------
        # Filtered QC
        # -----------------------------------------------------------------
        *(
            expand("outputs/{sample}/reports/QC/{sample}_NanoPlot-report.html", sample=SAMPLES)
            if PREPROCESSING_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/QC/{sample}_nanoQC.html", sample=SAMPLES)
            if PREPROCESSING_ENABLED else []
        ),

        # -----------------------------------------------------------------
        # Profiling branch
        # -----------------------------------------------------------------
        *(
            expand("outputs/{sample}/reports/taxonomy/kraken2/{sample}.kraken2.output.txt", sample=SAMPLES)
            if TAXONOMY_ENABLED and "kraken2" in TAX_TOOLS else []
        ),
        *(
            expand("outputs/{sample}/reports/taxonomy/kraken2/{sample}.kraken2.report.txt", sample=SAMPLES)
            if TAXONOMY_ENABLED and "kraken2" in TAX_TOOLS else []
        ),
        *(
            expand("outputs/{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.classification.tsv", sample=SAMPLES)
            if TAXONOMY_ENABLED and "centrifuger" in TAX_TOOLS else []
        ),
        *(
            expand("outputs/{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.report.tsv", sample=SAMPLES)
            if TAXONOMY_ENABLED and "centrifuger" in TAX_TOOLS else []
        ),
        *(
            expand("outputs/{sample}/stages/profile/{sample}.done", sample=SAMPLES)
            if RUN_PROFILE and RUN_ASSEMBLY and SERIAL_PROFILE_THEN_ASSEMBLY else []
        ),

        # -----------------------------------------------------------------
        # Assembly / binning / genome branch
        # -----------------------------------------------------------------
        *(
            expand("outputs/{sample}/reports/assembly/metaflye/{sample}/assembly.fasta", sample=SAMPLES)
            if ASSEMBLY_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/assembly/metaflye/{sample}/assembly_info.txt", sample=SAMPLES)
            if ASSEMBLY_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/assembly/metaflye/{sample}/assembly_graph.gfa", sample=SAMPLES)
            if ASSEMBLY_ENABLED else []
        ),

        # Assembly QC
        *(
            expand("outputs/{sample}/reports/assembly/assembly_qc/metaquast/{sample}/report.html", sample=SAMPLES)
            if ASSEMBLY_ENABLED and METAQUAST_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/assembly/assembly_qc/metaquast/{sample}/report.tsv", sample=SAMPLES)
            if ASSEMBLY_ENABLED and METAQUAST_ENABLED else []
        ),

        # Mapping and depth
        *(
            expand("outputs/{sample}/reports/mapping_depth/{sample}.vs_assembly.sorted.bam", sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/mapping_depth/{sample}.vs_assembly.sorted.bam.bai", sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/mapping_depth/{sample}.depth.txt", sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED else []
        ),

        # Binning
        *(
            expand("outputs/{sample}/reports/binning/metabat2/.done", sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "metabat2" in BINNING_TOOLS else []
        ),
        *(
            expand("outputs/{sample}/reports/binning/semibin2/.done", sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "semibin2" in BINNING_TOOLS else []
        ),
        *(
            expand("outputs/{sample}/reports/binning/comebin/.done", sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "comebin" in BINNING_TOOLS else []
        ),

        # Collected bins
        *(
            expand("outputs/{sample}/reports/binning/all_bins/.done", sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and len(BINNING_TOOLS) > 0 else []
        ),

        # Normalized contig-to-bin tables
        *(
            expand("outputs/{sample}/reports/binning/normalized/metabat2.contig2bin.tsv", sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "metabat2" in BINNING_TOOLS else []
        ),
        *(
            expand("outputs/{sample}/reports/binning/normalized/semibin2.contig2bin.tsv", sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "semibin2" in BINNING_TOOLS else []
        ),
        *(
            expand("outputs/{sample}/reports/binning/normalized/comebin.contig2bin.tsv", sample=SAMPLES)
            if ASSEMBLY_ENABLED and MAPPING_ENABLED and BINNING_ENABLED and "comebin" in BINNING_TOOLS else []
        ),

        # Prokaryotic refinement
        *(
            expand("outputs/{sample}/reports/refinement/prok/dastool/{sample}_DASTool_summary.tsv", sample=SAMPLES)
            if DASTOOL_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/prok/dastool/{sample}_DASTool_contig2bin.tsv", sample=SAMPLES)
            if DASTOOL_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/prok/dastool/{sample}_DASTool_bins", sample=SAMPLES)
            if DASTOOL_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/prok/dastool/.done", sample=SAMPLES)
            if DASTOOL_ENABLED else []
        ),

        # Prokaryotic quality estimation
        *(
            expand("outputs/{sample}/reports/refinement/prok/checkm2/quality_report.tsv", sample=SAMPLES)
            if CHECKM2_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/prok/checkm2/.done", sample=SAMPLES)
            if CHECKM2_ENABLED else []
        ),

        # Prokaryotic classification
        *(
            expand("outputs/{sample}/reports/refinement/prok/gtdbtk/gtdbtk.bac120.summary.tsv", sample=SAMPLES)
            if GTDBTK_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/prok/gtdbtk/gtdbtk.ar53.summary.tsv", sample=SAMPLES)
            if GTDBTK_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/prok/gtdbtk/.done", sample=SAMPLES)
            if GTDBTK_ENABLED else []
        ),

        # Eukaryotic refinement
        *(
            expand("outputs/{sample}/reports/refinement/euk/acr/.done", sample=SAMPLES)
            if ACR_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/euk/euk_bins/kept_bins.tsv", sample=SAMPLES)
            if ACR_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/euk/euk_bins/.done", sample=SAMPLES)
            if ACR_ENABLED else []
        ),

        # Eukaryotic dereplication
        *(
            expand("outputs/{sample}/reports/refinement/euk/drep/data_tables/Cdb.csv", sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/euk/drep/.done", sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED else []
        ),

        # Eukaryotic quality estimation
        *(
            expand("outputs/{sample}/reports/refinement/euk/eukcc/eukcc.csv", sample=SAMPLES)
            if ACR_ENABLED and EUKCC_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/euk/eukcc/.done", sample=SAMPLES)
            if ACR_ENABLED and EUKCC_ENABLED else []
        ),

        # Final eukaryotic bin selection
        *(
            expand("outputs/{sample}/reports/refinement/euk/final_bins/selected_bins.tsv", sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/euk/final_bins/.done", sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED else []
        ),

        # Eukaryotic bin classification
        *(
            expand("outputs/{sample}/reports/refinement/euk/bat/bin2classification.txt", sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED and BAT_ENABLED else []
        ),
        *(
            expand("outputs/{sample}/reports/refinement/euk/bat/.done", sample=SAMPLES)
            if ACR_ENABLED and DREP_ENABLED and EUKCC_ENABLED and BAT_ENABLED else []
        ),

        # Genome-level final reports
        *(
            expand("outputs/{sample}/reports/final/{sample}.genome_inventory.tsv", sample=SAMPLES)
            if RUN_ASSEMBLY else []
        ),
        *(
            expand("outputs/{sample}/reports/final/{sample}.genome_summary.tsv", sample=SAMPLES)
            if RUN_ASSEMBLY else []
        ),
        *(
            expand("outputs/{sample}/reports/final/{sample}.bin_trace.tsv", sample=SAMPLES)
            if RUN_ASSEMBLY else []
        ),

        # Combined pipeline summary
        *(
            expand("outputs/{sample}/reports/final/{sample}.pipeline_summary.tsv", sample=SAMPLES)
            if RUN_PROFILE and RUN_ASSEMBLY else []
        )
