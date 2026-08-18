# -------------------------------------------------------------------------
# DAS Tool refinement
# -------------------------------------------------------------------------
rule dastool_refine:
    container:
        CONTAINERS["prok"]
    input:
        assembly="{sample}/reports/assembly/metaflye/{sample}/assembly.fasta",
        bin_tables=dastool_inputs
    output:
        summary="{sample}/reports/refinement/prok/dastool/{sample}_DASTool_summary.tsv",
        contig2bin="{sample}/reports/refinement/prok/dastool/{sample}_DASTool_contig2bin.tsv",
        bins=directory("{sample}/reports/refinement/prok/dastool/{sample}_DASTool_bins"),
        done="{sample}/reports/refinement/prok/dastool/.done"
    threads:
        P["threads"].get("dastool", 16)
    params:
        outdir="{sample}/reports/refinement/prok/dastool",
        out_prefix="{sample}/reports/refinement/prok/dastool/{sample}",
        labels=lambda wc: dastool_labels(),
        input_string=lambda wc: dastool_input_string(wc),
        extra=DASTOOL_EXTRA,
        score_threshold=DASTOOL_SCORE_THRESHOLD,
        search_engine=DASTOOL_SEARCH_ENGINE,
        write_bin_evals="--write_bin_evals" if DASTOOL_WRITE_BIN_EVALS else "",
        write_unbinned="--write_unbinned" if DASTOOL_WRITE_UNBINNED else "",
        db_directory=DASTOOL_DB_DIRECTORY,
        mapping_enabled=MAPPING_ENABLED,
        binning_enabled=BINNING_ENABLED,
        dastool_enabled=DASTOOL_ENABLED,
        tools_count=len(enabled_dastool_tools())
    shell:
        r"""
        python - <<'PY'
import sys
mapping_enabled = "{params.mapping_enabled}"
binning_enabled = "{params.binning_enabled}"
dastool_enabled = "{params.dastool_enabled}"
tools_count = int("{params.tools_count}")

if mapping_enabled != "True":
    sys.stderr.write("ERROR: DAS Tool refinement requires mapping.enabled: true.\n")
    sys.exit(1)
if binning_enabled != "True":
    sys.stderr.write("ERROR: DAS Tool refinement requires binning.enabled: true.\n")
    sys.exit(1)
if dastool_enabled != "True":
    sys.stderr.write("ERROR: refinement.prok.dastool is false, but dastool_refine was requested.\n")
    sys.exit(1)
if tools_count == 0:
    sys.stderr.write("ERROR: DAS Tool refinement requires at least one enabled binner.\n")
    sys.exit(1)
PY
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        DAS_Tool \
          -i {params.input_string} \
          -l {params.labels} \
          -c {input.assembly} \
          -o {params.out_prefix} \
          --search_engine {params.search_engine} \
          --score_threshold {params.score_threshold} \
          --dbDirectory {params.db_directory} \
          --write_bins \
          {params.write_bin_evals} \
          {params.write_unbinned} \
          -t {threads} \
          {params.extra}

        touch {output.done}
        """

# -------------------------------------------------------------------------
# CheckM2 quality assessment
# -------------------------------------------------------------------------
rule checkm2_prok:
    container:
        CONTAINERS["prok"]
    input:
        bins="{sample}/reports/refinement/prok/dastool/{sample}_DASTool_bins",
        dastool_done="{sample}/reports/refinement/prok/dastool/.done"
    output:
        quality="{sample}/reports/refinement/prok/checkm2/quality_report.tsv",
        done="{sample}/reports/refinement/prok/checkm2/.done"
    threads:
        P["threads"].get("checkm2", 24)
    params:
        outdir="{sample}/reports/refinement/prok/checkm2",
        db=CHECKM2_DB,
        extension=CHECKM2_EXTENSION,
        extra=CHECKM2_EXTRA,
        mapping_enabled=MAPPING_ENABLED,
        binning_enabled=BINNING_ENABLED,
        dastool_enabled=DASTOOL_ENABLED,
        checkm2_enabled=CHECKM2_ENABLED
    shell:
        r"""
        python - <<'PY'
import sys
mapping_enabled = "{params.mapping_enabled}"
binning_enabled = "{params.binning_enabled}"
dastool_enabled = "{params.dastool_enabled}"
checkm2_enabled = "{params.checkm2_enabled}"

if mapping_enabled != "True":
    sys.stderr.write("ERROR: CheckM2 requires mapping.enabled: true.\n")
    sys.exit(1)
if binning_enabled != "True":
    sys.stderr.write("ERROR: CheckM2 requires binning.enabled: true.\n")
    sys.exit(1)
if dastool_enabled != "True":
    sys.stderr.write("ERROR: CheckM2 requires refinement.prok.dastool: true.\n")
    sys.exit(1)
if checkm2_enabled != "True":
    sys.stderr.write("ERROR: refinement.prok.checkm2 is false, but checkm2_prok was requested.\n")
    sys.exit(1)
PY
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        checkm2 predict \
          --input {input.bins} \
          --output-directory {params.outdir} \
          --threads {threads} \
          -x {params.extension} \
          --database_path {params.db} \
          {params.extra}

        touch {output.done}
        """

# -------------------------------------------------------------------------
# GTDB-Tk classification
# -------------------------------------------------------------------------
rule gtdbtk_classify:
    container:
        CONTAINERS["prok"]
    input:
        bins="{sample}/reports/refinement/prok/dastool/{sample}_DASTool_bins",
        checkm2_done="{sample}/reports/refinement/prok/checkm2/.done"
    output:
        bac120="{sample}/reports/refinement/prok/gtdbtk/gtdbtk.bac120.summary.tsv",
        ar53="{sample}/reports/refinement/prok/gtdbtk/gtdbtk.ar53.summary.tsv",
        done="{sample}/reports/refinement/prok/gtdbtk/.done"
    threads:
        P["threads"].get("gtdbtk", 32)
    params:
        outdir="{sample}/reports/refinement/prok/gtdbtk",
        db_env=(lambda wc: f'export GTDBTK_DATA_PATH="{GTDBTK_DB}";' if GTDBTK_DB else ""),
        mash_arg=(lambda wc: f"--mash_db {GTDBTK_MASH_DB}" if GTDBTK_MASH_DB else ""),
        extra=GTDBTK_EXTRA,
        mapping_enabled=MAPPING_ENABLED,
        binning_enabled=BINNING_ENABLED,
        dastool_enabled=DASTOOL_ENABLED,
        checkm2_enabled=CHECKM2_ENABLED,
        gtdbtk_enabled=GTDBTK_ENABLED
    shell:
        r"""
        python - <<'PY'
import sys
mapping_enabled = "{params.mapping_enabled}"
binning_enabled = "{params.binning_enabled}"
dastool_enabled = "{params.dastool_enabled}"
checkm2_enabled = "{params.checkm2_enabled}"
gtdbtk_enabled = "{params.gtdbtk_enabled}"

if mapping_enabled != "True":
    sys.stderr.write("ERROR: GTDB-Tk requires mapping.enabled: true.\n")
    sys.exit(1)
if binning_enabled != "True":
    sys.stderr.write("ERROR: GTDB-Tk requires binning.enabled: true.\n")
    sys.exit(1)
if dastool_enabled != "True":
    sys.stderr.write("ERROR: GTDB-Tk requires refinement.prok.dastool: true.\n")
    sys.exit(1)
if checkm2_enabled != "True":
    sys.stderr.write("ERROR: GTDB-Tk requires refinement.prok.checkm2: true.\n")
    sys.exit(1)
if gtdbtk_enabled != "True":
    sys.stderr.write("ERROR: refinement.prok.gtdbtk is false, but gtdbtk_classify was requested.\n")
    sys.exit(1)
PY
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        {params.db_env}
        gtdbtk classify_wf \
          --genome_dir {input.bins} \
          --out_dir {params.outdir} \
          --cpus {threads} \
          --extension fa \
          {params.mash_arg} \
          {params.extra}

        BAC_SRC=""
        AR_SRC=""

        if [ -f "{params.outdir}/gtdbtk.bac120.summary.tsv" ]; then
            BAC_SRC="{params.outdir}/gtdbtk.bac120.summary.tsv"
        elif [ -f "{params.outdir}/classify/gtdbtk.bac120.summary.tsv" ]; then
            BAC_SRC="{params.outdir}/classify/gtdbtk.bac120.summary.tsv"
        elif [ -f "{params.outdir}/classify/ani_screen/gtdbtk.bac120.ani_summary.tsv" ]; then
            BAC_SRC="{params.outdir}/classify/ani_screen/gtdbtk.bac120.ani_summary.tsv"
        fi

        if [ -f "{params.outdir}/gtdbtk.ar53.summary.tsv" ]; then
            AR_SRC="{params.outdir}/gtdbtk.ar53.summary.tsv"
        elif [ -f "{params.outdir}/classify/gtdbtk.ar53.summary.tsv" ]; then
            AR_SRC="{params.outdir}/classify/gtdbtk.ar53.summary.tsv"
        elif [ -f "{params.outdir}/classify/ani_screen/gtdbtk.ar53.ani_summary.tsv" ]; then
            AR_SRC="{params.outdir}/classify/ani_screen/gtdbtk.ar53.ani_summary.tsv"
        fi

        if [ -n "$BAC_SRC" ]; then
            if [ "$(readlink -f "$BAC_SRC")" != "$(readlink -f "{output.bac120}")" ]; then
                cp "$BAC_SRC" "{output.bac120}"
            fi
        else
            printf "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings\n" > "{output.bac120}"
        fi

        if [ -n "$AR_SRC" ]; then
            if [ "$(readlink -f "$AR_SRC")" != "$(readlink -f "{output.ar53}")" ]; then
                cp "$AR_SRC" "{output.ar53}"
            fi
        else
            printf "user_genome\tclassification\tfastani_reference\tfastani_reference_radius\tfastani_taxonomy\tfastani_ani\tfastani_af\tclosest_placement_reference\tclosest_placement_radius\tclosest_placement_taxonomy\tclosest_placement_ani\tclosest_placement_af\tpplacer_taxonomy\tclassification_method\tnote\tother_related_references(genome_id,species_name,radius,ANI,AF)\tmsa_percent\ttranslation_table\tred_value\twarnings\n" > "{output.ar53}"
        fi

        touch {output.done}
        """

