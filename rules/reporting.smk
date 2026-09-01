# -------------------------------------------------------------------------
# Versions manifest in the report folder
# Written by the run wrapper before Snakemake starts.
# This rule serves as a fallback when Snakemake is called directly.
# -------------------------------------------------------------------------
_PIPELINE_ROOT = os.path.dirname(config.get("path_scripts", os.path.join(PIPELINE_DIR, "scripts")))

rule write_versions_manifest:
    output:
        os.path.join("Reports", "versions.json")
    params:
        src=os.path.join(_PIPELINE_ROOT, "metadata", "versions.json")
    shell:
        r"""
        cp {params.src} {output}
        """


# -------------------------------------------------------------------------
# Parse NanoPlot metrics → JSON (raw + filtered)
# -------------------------------------------------------------------------
rule parse_nanoplot_metrics:
    input:
        raw_stats=qc_path("Raw", "{sample}", "{sample}_NanoStats.txt"),
        filtered_stats=qc_path("Filtered", "{sample}", "{sample}_NanoStats.txt")
    output:
        json=report_path("{sample}", "{sample}.nanoplot_metrics.json")
    shell:
        r"""
        mkdir -p "$(dirname {output.json})"
        python {SCRIPTS_DIR}/parse_nanoplot_metrics.py \
          --raw-nanostats "{input.raw_stats}" \
          --filtered-nanostats "{input.filtered_stats}" \
          --output {output.json}
        """


# -------------------------------------------------------------------------
# Parse taxonomy reports → JSON (top-N species per classifier)
# -------------------------------------------------------------------------
rule parse_taxonomy_report:
    input:
        kraken2_report=lambda wc: (
            taxonomy_path("Kraken2", wc.sample, f"{wc.sample}.kraken2.report.txt")
            if TAXONOMY_ENABLED and "kraken2" in TAX_TOOLS else []
        ),
        centrifuger_report=lambda wc: (
            taxonomy_path("Centrifuger", wc.sample, f"{wc.sample}.centrifuger.report.tsv")
            if TAXONOMY_ENABLED and "centrifuger" in TAX_TOOLS else []
        )
    output:
        json=report_path("{sample}", "{sample}.taxonomy_top10.json")
    shell:
        r"""
        mkdir -p "$(dirname {output.json})"
        python {SCRIPTS_DIR}/parse_taxonomy.py \
          --kraken2-report "{input.kraken2_report}" \
          --centrifuger-report "{input.centrifuger_report}" \
          --top-n 10 \
          --output {output.json}
        """


# -------------------------------------------------------------------------
# Parse assembly QC metrics (N50, contig count, ...) → JSON
#
# Assembly runs for every sample in an assembly-enabled run mode, regardless
# of whether it later passes the binning gate, so this is always computed
# per-sample when RUN_ASSEMBLY is true (unlike genome/pipeline summaries,
# which only exist for samples that were binned).
# -------------------------------------------------------------------------
rule parse_assembly_stats:
    input:
        assembly_info=assembly_path("MetaFlye", "{sample}", "assembly_info.txt"),
        metaquast_tsv=(
            assembly_path("MetaQUAST", "{sample}", "report.tsv")
            if METAQUAST_ENABLED else []
        )
    output:
        json=report_path("{sample}", "{sample}.assembly_stats.json")
    params:
        metaquast_tsv=assembly_path("MetaQUAST", "{sample}", "report.tsv") if METAQUAST_ENABLED else ""
    shell:
        r"""
        mkdir -p "$(dirname {output.json})"
        python {SCRIPTS_DIR}/parse_assembly_stats.py \
          --sample {wildcards.sample} \
          --assembly-info "{input.assembly_info}" \
          --metaquast-report "{params.metaquast_tsv}" \
          --output {output.json}
        """


# -------------------------------------------------------------------------
# Render HTML report (single run-level report for all samples)
# -------------------------------------------------------------------------
rule render_run_report:
    input:
        taxonomy_jsons=expand(
            report_path("{sample}", "{sample}.taxonomy_top10.json"),
            sample=SAMPLES
        ) if TAXONOMY_ENABLED else [],
        nanoplot_jsons=expand(
            report_path("{sample}", "{sample}.nanoplot_metrics.json"),
            sample=SAMPLES
        ),
        assembly_stats_jsons=expand(
            report_path("{sample}", "{sample}.assembly_stats.json"),
            sample=SAMPLES
        ) if RUN_ASSEMBLY else [],
        # Only samples that passed the post-assembly binning gate have
        # genome/pipeline summaries (see binning_gate_passed() in common.smk).
        # Each row/record embeds its own "sample" field, so render_report.py
        # matches them back to samples by name rather than by list position.
        genome_summaries=(
            lambda wc: [report_path(s, f"{s}.genome_summary.tsv") for s in binning_samples()]
        ) if RUN_ASSEMBLY else [],
        gtdbtk_bac120s=(
            lambda wc: [
                refinement_path("Prokaryotic", "GTDBTk", s, "gtdbtk.bac120.summary.tsv")
                for s in binning_samples()
            ]
        ) if GTDBTK_ENABLED else [],
        gtdbtk_ar53s=(
            lambda wc: [
                refinement_path("Prokaryotic", "GTDBTk", s, "gtdbtk.ar53.summary.tsv")
                for s in binning_samples()
            ]
        ) if GTDBTK_ENABLED else [],
        genome_inventories=(
            lambda wc: [report_path(s, f"{s}.genome_inventory.tsv") for s in binning_samples()]
        ) if RUN_ASSEMBLY else [],
        pipeline_summaries=(
            lambda wc: [report_path(s, f"{s}.pipeline_summary.tsv") for s in binning_samples()]
        ) if RUN_PROFILE and RUN_ASSEMBLY else [],
        host_stats=expand(
            preprocessing_path("{sample}", "{sample}.host_removal_stats.json"),
            sample=SAMPLES
        ) if HOST_REMOVAL_ENABLED else [],
        versions_json=os.path.join("Reports", "versions.json")
    output:
        html=os.path.join("Reports", "report.html")
    params:
        samples=",".join(SAMPLES),
        binning_enabled=BINNING_ENABLED,
        run_mode=RUN_MODE,
        filtering_method=FILTERING_METHOD,
        preprocessing_enabled=PREPROCESSING_ENABLED,
        host_removal_enabled=HOST_REMOVAL_ENABLED,
        host_ref=HOST_REF,
        chopper_min_length=P["chopper"]["min_length"],
        chopper_quality=P["chopper"]["quality_threshold"],
        filtlong_min_length=P["filtlong"]["min_length"],
        filtlong_keep_percent=P["filtlong"]["keep_percent"],
        workflow_png=os.path.join(
            os.path.dirname(config.get("path_scripts", "scripts")), "_workflow_.png"
        ),
        plotly_js=os.path.join(
            os.path.dirname(config.get("path_scripts", "scripts")),
            "templates", "report", "plotly-v1.58.5.js"
        )
    shell:
        r"""
        python {SCRIPTS_DIR}/render_report.py \
          --samples "{params.samples}" \
          --run-mode {params.run_mode} \
          --binning-enabled "{params.binning_enabled}" \
          --taxonomy-jsons "{input.taxonomy_jsons}" \
          --nanoplot-jsons "{input.nanoplot_jsons}" \
          --assembly-stats-jsons "{input.assembly_stats_jsons}" \
          --genome-summaries "{input.genome_summaries}" \
          --gtdbtk-bac120s "{input.gtdbtk_bac120s}" \
          --gtdbtk-ar53s "{input.gtdbtk_ar53s}" \
          --genome-inventories "{input.genome_inventories}" \
          --pipeline-summaries "{input.pipeline_summaries}" \
          --host-stats-jsons "{input.host_stats}" \
          --versions-json "{input.versions_json}" \
          --filtering-method "{params.filtering_method}" \
          --preprocessing-enabled "{params.preprocessing_enabled}" \
          --host-removal-enabled "{params.host_removal_enabled}" \
          --host-ref "{params.host_ref}" \
          --chopper-min-length "{params.chopper_min_length}" \
          --chopper-quality "{params.chopper_quality}" \
          --filtlong-min-length "{params.filtlong_min_length}" \
          --filtlong-keep-percent "{params.filtlong_keep_percent}" \
          --workflow-png "{params.workflow_png}" \
          --plotly-js "{params.plotly_js}" \
          --output {output.html}
        """

# -------------------------------------------------------------------------
# Branch 2 genome reporting
# -------------------------------------------------------------------------
rule report_branch2_genomes:
    input:
        prok_bins=refinement_path("Prokaryotic", "DASTool", "{sample}", "{sample}_DASTool_bins"),
        checkm2=refinement_path("Prokaryotic", "CheckM2", "{sample}", "quality_report.tsv"),
        gtdb_bac=refinement_path("Prokaryotic", "GTDBTk", "{sample}", "gtdbtk.bac120.summary.tsv"),
        gtdb_ar=refinement_path("Prokaryotic", "GTDBTk", "{sample}", "gtdbtk.ar53.summary.tsv"),
        dastool_map=refinement_path("Prokaryotic", "DASTool", "{sample}", "{sample}_DASTool_contig2bin.tsv"),
        euk_bins=refinement_path("Eukaryotic", "FinalBins", "{sample}", "selected_bins"),
        eukcc=refinement_path("Eukaryotic", "EukCC", "{sample}", "eukcc.csv"),
        bat=refinement_path("Eukaryotic", "BAT", "{sample}", "bin2classification.txt"),
        drep_clusters=refinement_path("Eukaryotic", "dRep", "{sample}", "data_tables", "Cdb.csv"),
        kept_manifest=refinement_path("Eukaryotic", "EukBins", "{sample}", "kept_bins.tsv"),
        selected_manifest=refinement_path("Eukaryotic", "FinalBins", "{sample}", "selected_bins.tsv")
    output:
        inventory=report_path("{sample}", "{sample}.genome_inventory.tsv"),
        summary=report_path("{sample}", "{sample}.genome_summary.tsv"),
        trace=report_path("{sample}", "{sample}.bin_trace.tsv")
    shell:
        r"""
        mkdir -p "$(dirname {output.inventory})"
        python {SCRIPTS_DIR}/merge_final_reports.py \
          --sample {wildcards.sample} \
          --prok-bins {input.prok_bins} \
          --checkm2 {input.checkm2} \
          --gtdb-bac {input.gtdb_bac} \
          --gtdb-ar {input.gtdb_ar} \
          --dastool-map {input.dastool_map} \
          --euk-bins {input.euk_bins} \
          --eukcc {input.eukcc} \
          --bat {input.bat} \
          --drep-clusters {input.drep_clusters} \
          --kept-manifest {input.kept_manifest} \
          --selected-manifest {input.selected_manifest} \
          --inventory {output.inventory} \
          --summary {output.summary} \
          --trace {output.trace}
        """

# -------------------------------------------------------------------------
# Profiling stage barrier
# -------------------------------------------------------------------------
rule profile_stage_complete:
    input:
        kraken2_out=taxonomy_path("Kraken2", "{sample}", "{sample}.kraken2.output.txt") if "kraken2" in TAX_TOOLS else [],
        kraken2_report=taxonomy_path("Kraken2", "{sample}", "{sample}.kraken2.report.txt") if "kraken2" in TAX_TOOLS else [],
        centrifuger_classif=taxonomy_path("Centrifuger", "{sample}", "{sample}.centrifuger.classification.tsv") if "centrifuger" in TAX_TOOLS else [],
        centrifuger_report=taxonomy_path("Centrifuger", "{sample}", "{sample}.centrifuger.report.tsv") if "centrifuger" in TAX_TOOLS else []
    output:
        done=stage_path("Profile", "{sample}", "{sample}.done")
    shell:
        r"""
        mkdir -p "$(dirname {output.done})"
        touch {output.done}
        """

# -------------------------------------------------------------------------
# Top-level pipeline summary
# -------------------------------------------------------------------------
rule report_pipeline_summary:
    input:
        kraken2_report=taxonomy_path("Kraken2", "{sample}", "{sample}.kraken2.report.txt") if "kraken2" in TAX_TOOLS else [],
        centrifuger_report=taxonomy_path("Centrifuger", "{sample}", "{sample}.centrifuger.report.tsv") if "centrifuger" in TAX_TOOLS else [],
        genome_summary=report_path("{sample}", "{sample}.genome_summary.tsv"),
        inventory=report_path("{sample}", "{sample}.genome_inventory.tsv")
    output:
        summary=report_path("{sample}", "{sample}.pipeline_summary.tsv")
    shell:
        r"""
        python {SCRIPTS_DIR}/report_pipeline_summary.py \
          --sample {wildcards.sample} \
          --kraken2-report "{input.kraken2_report}" \
          --centrifuger-report "{input.centrifuger_report}" \
          --genome-summary {input.genome_summary} \
          --inventory {input.inventory} \
          --output {output.summary}
        """
