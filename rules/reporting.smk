# -------------------------------------------------------------------------
# Versions manifest at output root
# Written by the run wrapper before Snakemake starts.
# This rule serves as a fallback when Snakemake is called directly.
# -------------------------------------------------------------------------
_PIPELINE_ROOT = os.path.dirname(config.get("path_scripts", os.path.join(PIPELINE_DIR, "scripts")))

rule write_versions_manifest:
    output:
        "versions.json"
    params:
        src=os.path.join(_PIPELINE_ROOT, "metadata", "versions.json")
    shell:
        r"""
        cp {params.src} {output}
        """


# -------------------------------------------------------------------------
# Parse taxonomy reports → JSON (top-N species per classifier)
# -------------------------------------------------------------------------
rule parse_taxonomy_report:
    input:
        kraken2_report=lambda wc: (
            f"{wc.sample}/reports/taxonomy/kraken2/{wc.sample}.kraken2.report.txt"
            if TAXONOMY_ENABLED and "kraken2" in TAX_TOOLS else []
        ),
        centrifuger_report=lambda wc: (
            f"{wc.sample}/reports/taxonomy/centrifuger/{wc.sample}.centrifuger.report.tsv"
            if TAXONOMY_ENABLED and "centrifuger" in TAX_TOOLS else []
        )
    output:
        json="{sample}/reports/final/{sample}.taxonomy_top10.json"
    shell:
        r"""
        mkdir -p {wildcards.sample}/reports/final
        python {SCRIPTS_DIR}/parse_taxonomy.py \
          --kraken2-report "{input.kraken2_report}" \
          --centrifuger-report "{input.centrifuger_report}" \
          --top-n 10 \
          --output {output.json}
        """


# -------------------------------------------------------------------------
# Render HTML report (single run-level report for all samples)
# -------------------------------------------------------------------------
rule render_run_report:
    input:
        taxonomy_jsons=expand(
            "{sample}/reports/final/{sample}.taxonomy_top10.json",
            sample=SAMPLES
        ) if TAXONOMY_ENABLED else [],
        genome_summaries=expand(
            "{sample}/reports/final/{sample}.genome_summary.tsv",
            sample=SAMPLES
        ) if RUN_ASSEMBLY else [],
        genome_inventories=expand(
            "{sample}/reports/final/{sample}.genome_inventory.tsv",
            sample=SAMPLES
        ) if RUN_ASSEMBLY else [],
        pipeline_summaries=expand(
            "{sample}/reports/final/{sample}.pipeline_summary.tsv",
            sample=SAMPLES
        ) if RUN_PROFILE and RUN_ASSEMBLY else [],
        host_stats=expand(
            "{sample}/reports/preprocessing/{sample}.host_removal_stats.json",
            sample=SAMPLES
        ) if HOST_REMOVAL_ENABLED else [],
        versions_json="versions.json"
    output:
        html="report.html"
    params:
        samples=",".join(SAMPLES),
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
        )
    shell:
        r"""
        python {SCRIPTS_DIR}/render_report.py \
          --samples "{params.samples}" \
          --run-mode {params.run_mode} \
          --taxonomy-jsons "{input.taxonomy_jsons}" \
          --genome-summaries "{input.genome_summaries}" \
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
          --output {output.html}
        """

# -------------------------------------------------------------------------
# Branch 2 genome reporting
# -------------------------------------------------------------------------
rule report_branch2_genomes:
    input:
        prok_bins="{sample}/reports/refinement/prok/dastool/{sample}_DASTool_bins",
        checkm2="{sample}/reports/refinement/prok/checkm2/quality_report.tsv",
        gtdb_bac="{sample}/reports/refinement/prok/gtdbtk/gtdbtk.bac120.summary.tsv",
        gtdb_ar="{sample}/reports/refinement/prok/gtdbtk/gtdbtk.ar53.summary.tsv",
        dastool_map="{sample}/reports/refinement/prok/dastool/{sample}_DASTool_contig2bin.tsv",
        euk_bins="{sample}/reports/refinement/euk/final_bins/selected_bins",
        eukcc="{sample}/reports/refinement/euk/eukcc/eukcc.csv",
        bat="{sample}/reports/refinement/euk/bat/bin2classification.txt",
        drep_clusters="{sample}/reports/refinement/euk/drep/data_tables/Cdb.csv",
        kept_manifest="{sample}/reports/refinement/euk/euk_bins/kept_bins.tsv",
        selected_manifest="{sample}/reports/refinement/euk/final_bins/selected_bins.tsv"
    output:
        inventory="{sample}/reports/final/{sample}.genome_inventory.tsv",
        summary="{sample}/reports/final/{sample}.genome_summary.tsv",
        trace="{sample}/reports/final/{sample}.bin_trace.tsv"
    shell:
        r"""
        mkdir -p {wildcards.sample}/reports/final
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
        kraken2_out="{sample}/reports/taxonomy/kraken2/{sample}.kraken2.output.txt" if "kraken2" in TAX_TOOLS else [],
        kraken2_report="{sample}/reports/taxonomy/kraken2/{sample}.kraken2.report.txt" if "kraken2" in TAX_TOOLS else [],
        centrifuger_classif="{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.classification.tsv" if "centrifuger" in TAX_TOOLS else [],
        centrifuger_report="{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.report.tsv" if "centrifuger" in TAX_TOOLS else []
    output:
        done="{sample}/stages/profile/{sample}.done"
    shell:
        r"""
        mkdir -p {wildcards.sample}/stages/profile
        touch {output.done}
        """

# -------------------------------------------------------------------------
# Top-level pipeline summary
# -------------------------------------------------------------------------
rule report_pipeline_summary:
    input:
        kraken2_report="{sample}/reports/taxonomy/kraken2/{sample}.kraken2.report.txt" if "kraken2" in TAX_TOOLS else [],
        centrifuger_report="{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.report.tsv" if "centrifuger" in TAX_TOOLS else [],
        genome_summary="{sample}/reports/final/{sample}.genome_summary.tsv",
        inventory="{sample}/reports/final/{sample}.genome_inventory.tsv"
    output:
        summary="{sample}/reports/final/{sample}.pipeline_summary.tsv"
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
