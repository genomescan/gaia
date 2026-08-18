
# -------------------------------------------------------------------------
# Parse taxonomy reports → JSON (top-N species per classifier)
# -------------------------------------------------------------------------
rule parse_taxonomy_report:
    input:
        kraken2_report=lambda wc: (
            f"outputs/{wc.sample}/reports/taxonomy/kraken2/{wc.sample}.kraken2.report.txt"
            if TAXONOMY_ENABLED and "kraken2" in TAX_TOOLS else []
        ),
        centrifuger_report=lambda wc: (
            f"outputs/{wc.sample}/reports/taxonomy/centrifuger/{wc.sample}.centrifuger.report.tsv"
            if TAXONOMY_ENABLED and "centrifuger" in TAX_TOOLS else []
        )
    output:
        json="outputs/{sample}/reports/final/{sample}.taxonomy_top10.json"
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/final
        python {SCRIPTS_DIR}/parse_taxonomy.py \
          --kraken2-report "{input.kraken2_report}" \
          --centrifuger-report "{input.centrifuger_report}" \
          --top-n 10 \
          --output {output.json}
        """

# -------------------------------------------------------------------------
# Render HTML report
# -------------------------------------------------------------------------
rule render_report:
    input:
        taxonomy_json=lambda wc: (
            f"outputs/{wc.sample}/reports/final/{wc.sample}.taxonomy_top10.json"
            if TAXONOMY_ENABLED else []
        ),
        genome_summary=lambda wc: (
            f"outputs/{wc.sample}/reports/final/{wc.sample}.genome_summary.tsv"
            if RUN_ASSEMBLY else []
        ),
        genome_inventory=lambda wc: (
            f"outputs/{wc.sample}/reports/final/{wc.sample}.genome_inventory.tsv"
            if RUN_ASSEMBLY else []
        ),
        pipeline_summary=lambda wc: (
            f"outputs/{wc.sample}/reports/final/{wc.sample}.pipeline_summary.tsv"
            if RUN_PROFILE and RUN_ASSEMBLY else []
        )
    output:
        html="outputs/{sample}/reports/final/{sample}.report.html"
    params:
        run_mode=RUN_MODE
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/final
        python {SCRIPTS_DIR}/render_report.py \
          --sample {wildcards.sample} \
          --run-mode {params.run_mode} \
          --taxonomy-json "{input.taxonomy_json}" \
          --genome-summary "{input.genome_summary}" \
          --genome-inventory "{input.genome_inventory}" \
          --pipeline-summary "{input.pipeline_summary}" \
          --output {output.html}
        """

# -------------------------------------------------------------------------
# Branch 2 genome reporting
# -------------------------------------------------------------------------
rule report_branch2_genomes:
    input:
        prok_bins="outputs/{sample}/reports/refinement/prok/dastool/{sample}_DASTool_bins",
        checkm2="outputs/{sample}/reports/refinement/prok/checkm2/quality_report.tsv",
        gtdb_bac="outputs/{sample}/reports/refinement/prok/gtdbtk/gtdbtk.bac120.summary.tsv",
        gtdb_ar="outputs/{sample}/reports/refinement/prok/gtdbtk/gtdbtk.ar53.summary.tsv",
        dastool_map="outputs/{sample}/reports/refinement/prok/dastool/{sample}_DASTool_contig2bin.tsv",
        euk_bins="outputs/{sample}/reports/refinement/euk/final_bins/selected_bins",
        eukcc="outputs/{sample}/reports/refinement/euk/eukcc/eukcc.csv",
        bat="outputs/{sample}/reports/refinement/euk/bat/bin2classification.txt",
        drep_clusters="outputs/{sample}/reports/refinement/euk/drep/data_tables/Cdb.csv",
        kept_manifest="outputs/{sample}/reports/refinement/euk/euk_bins/kept_bins.tsv",
        selected_manifest="outputs/{sample}/reports/refinement/euk/final_bins/selected_bins.tsv"
    output:
        inventory="outputs/{sample}/reports/final/{sample}.genome_inventory.tsv",
        summary="outputs/{sample}/reports/final/{sample}.genome_summary.tsv",
        trace="outputs/{sample}/reports/final/{sample}.bin_trace.tsv"
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/final
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
        kraken2_out="outputs/{sample}/reports/taxonomy/kraken2/{sample}.kraken2.output.txt" if "kraken2" in TAX_TOOLS else [],
        kraken2_report="outputs/{sample}/reports/taxonomy/kraken2/{sample}.kraken2.report.txt" if "kraken2" in TAX_TOOLS else [],
        centrifuger_classif="outputs/{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.classification.tsv" if "centrifuger" in TAX_TOOLS else [],
        centrifuger_report="outputs/{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.report.tsv" if "centrifuger" in TAX_TOOLS else []
    output:
        done="outputs/{sample}/stages/profile/{sample}.done"
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/stages/profile
        touch {output.done}
        """

# -------------------------------------------------------------------------
# Top-level pipeline summary
# -------------------------------------------------------------------------
rule report_pipeline_summary:
    input:
        kraken2_report="outputs/{sample}/reports/taxonomy/kraken2/{sample}.kraken2.report.txt" if "kraken2" in TAX_TOOLS else [],
        centrifuger_report="outputs/{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.report.tsv" if "centrifuger" in TAX_TOOLS else [],
        genome_summary="outputs/{sample}/reports/final/{sample}.genome_summary.tsv",
        inventory="outputs/{sample}/reports/final/{sample}.genome_inventory.tsv"
    output:
        summary="outputs/{sample}/reports/final/{sample}.pipeline_summary.tsv"
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
