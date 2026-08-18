rule acr_refine:
    container:
        CONTAINERS["euk"]
    input:
        bins_dir="{sample}/reports/binning/all_bins",
        coverage="{sample}/reports/mapping_depth/{sample}.depth.txt",
        bins_done="{sample}/reports/binning/all_bins/.done",
        acr_db_dir=config["refinement"]["euk"]["acr_params"]["db_path"]
    output:
        refined_dir=directory("{sample}/reports/refinement/euk/acr/refined_bins"),
        done="{sample}/reports/refinement/euk/acr/.done"
    threads:
        P["threads"].get("acr", 16)
    params:
        outdir="{sample}/reports/refinement/euk/acr",
        prefix=ACR_PREFIX,
        min_size=ACR_MIN_SIZE,
        bypass=ACR_BYPASS,
        from_jgi=ACR_FROM_JGI,
        run_gmes=ACR_RUN_GMESEUK,
        target=ACR_TARGET,
        comp=ACR_COMP,
        cont=ACR_CONT,
        extra=ACR_EXTRA,
        mapping_enabled=MAPPING_ENABLED,
        binning_enabled=BINNING_ENABLED,
        acr_enabled=ACR_ENABLED
    shell:
        r"""
        python - <<'PY'
import sys
import pathlib

mapping_enabled = "{params.mapping_enabled}"
binning_enabled = "{params.binning_enabled}"
acr_enabled = "{params.acr_enabled}"

host_acr_db = pathlib.Path(r"{input.acr_db_dir}")

if mapping_enabled != "True":
    sys.stderr.write("ERROR: ACR requires mapping.enabled: true.\n")
    sys.exit(1)
if binning_enabled != "True":
    sys.stderr.write("ERROR: ACR requires binning.enabled: true.\n")
    sys.exit(1)
if acr_enabled != "True":
    sys.stderr.write("ERROR: refinement.euk.acr is false, but acr_refine was requested.\n")
    sys.exit(1)

if not host_acr_db.exists():
    sys.stderr.write(f"ERROR: Host ACR database directory does not exist: {{host_acr_db}}\n")
    sys.exit(1)

if not (host_acr_db / "eukcc").exists():
    sys.stderr.write(f"ERROR: Missing ACR eukcc database directory: {{host_acr_db / 'eukcc'}}\n")
    sys.exit(1)

if not (host_acr_db / "gtdbtk86").exists():
    sys.stderr.write(f"ERROR: Missing ACR gtdbtk86 database directory: {{host_acr_db / 'gtdbtk86'}}\n")
    sys.exit(1)
PY

        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        mkdir -p {output.refined_dir}

        BINS_DIR="$(realpath {input.bins_dir})"
        COVERAGE_FILE="$(realpath {input.coverage})"
        OUTDIR="$(realpath {params.outdir})"
        REFINED_DIR="$(realpath {output.refined_dir})"
        ACR_DB_PATH="$(realpath {input.acr_db_dir})"

        export ACR_DB_PATH
        export EUKCC2_DB="${{ACR_DB_PATH}}/eukcc"
        export GTDBTK_DATA_PATH="${{ACR_DB_PATH}}/gtdbtk86"

        acr \
          -g "$BINS_DIR" \
          -c "$COVERAGE_FILE" \
          -o "$OUTDIR" \
          -e fa \
          -p {params.prefix} \
          -t {threads} \
          -b {params.bypass} \
          -j {params.from_jgi} \
          -m {params.run_gmes} \
          --target {params.target} \
          --comp {params.comp} \
          --cont {params.cont} \
          {params.extra}

        if [ -f "$OUTDIR/Bin_Stat.csv" ]; then
          cp "$OUTDIR"/{params.prefix}*.fa "$REFINED_DIR"/ 2>/dev/null || true
        fi

        n_bins=$(find "$REFINED_DIR" -maxdepth 1 -type f -name "*.fa" | wc -l)
        if [ "$n_bins" -eq 0 ]; then
          echo "WARNING: ACR produced no refined bins for this sample." >&2
        fi

        touch {output.done}
        """

# -------------------------------------------------------------------------
# Keep only ACR-labeled eukaryotic bins
# -------------------------------------------------------------------------
rule keep_euk_from_acr:
    input:
        refined_dir="{sample}/reports/refinement/euk/acr/refined_bins",
        acr_done="{sample}/reports/refinement/euk/acr/.done"
    output:
        euk_dir=directory("{sample}/reports/refinement/euk/euk_bins"),
        manifest="{sample}/reports/refinement/euk/euk_bins/kept_bins.tsv",
        done="{sample}/reports/refinement/euk/euk_bins/.done"
    shell:
        r"""
        rm -rf {output.euk_dir}
        mkdir -p {output.euk_dir}

        python {SCRIPTS_DIR}/keep_euk_bins.py \
          --input-dir {input.refined_dir} \
          --output-dir {output.euk_dir} \
          --manifest {output.manifest}

        touch {output.done}
        """

# -------------------------------------------------------------------------
# dRep clustering only
# -------------------------------------------------------------------------
rule drep_euk_compare:
    container:
        CONTAINERS["euk"]
    input:
        euk_done="{sample}/reports/refinement/euk/euk_bins/.done",
        euk_dir="{sample}/reports/refinement/euk/euk_bins"
    output:
        cluster_info="{sample}/reports/refinement/euk/drep/data_tables/Cdb.csv",
        done="{sample}/reports/refinement/euk/drep/.done"
    threads:
        P["threads"].get("drep", 16)
    params:
        outdir="{sample}/reports/refinement/euk/drep",
        extra=DREP_EXTRA,
        drep_enabled=DREP_ENABLED
    shell:
        r"""
        python - <<'PY'
import sys, shutil, pathlib
drep_enabled = "{params.drep_enabled}"

if drep_enabled != "True":
    sys.stderr.write("ERROR: refinement.euk.drep is false, but drep_euk_compare was requested.\n")
    sys.exit(1)

indir = pathlib.Path(r"{input.euk_dir}")
outdir = pathlib.Path(r"{params.outdir}")
staged = outdir / "input_genomes"
tables = outdir / "data_tables"
done = pathlib.Path(r"{output.done}")
cluster = pathlib.Path(r"{output.cluster_info}")
genome_list = outdir / "genomes_to_compare.txt"

if outdir.exists():
    shutil.rmtree(outdir)

staged.mkdir(parents=True, exist_ok=True)
tables.mkdir(parents=True, exist_ok=True)

fasta_exts = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
files = []
if indir.exists():
    for p in sorted(indir.iterdir()):
        name = p.name.lower()
        if p.is_file() and any(name.endswith(ext) for ext in fasta_exts):
            dest = staged / p.name
            shutil.copy2(p, dest)
            files.append(dest)

if len(files) == 0:
    with open(cluster, "w") as out:
        out.write("genome\tsecondary_cluster\n")
    done.touch()
    sys.exit(0)

with open(genome_list, "w") as out:
    for f in files:
        out.write(str(f) + "\n")
PY

        if [ ! -s {params.outdir}/genomes_to_compare.txt ]; then
            touch {output.done}
            exit 0
        fi

        dRep compare {params.outdir} \
          -g $(cat {params.outdir}/genomes_to_compare.txt) \
          -p {threads} \
          {params.extra}

        touch {output.done}
        """

# -------------------------------------------------------------------------
# EukCC on all candidate euk bins
# -------------------------------------------------------------------------
rule eukcc_euk_candidates:
    container:
        CONTAINERS["euk"]
    input:
        euk_done="{sample}/reports/refinement/euk/euk_bins/.done",
        euk_dir="{sample}/reports/refinement/euk/euk_bins"
    output:
        summary="{sample}/reports/refinement/euk/eukcc/eukcc.csv",
        done="{sample}/reports/refinement/euk/eukcc/.done"
    threads:
        P["threads"].get("eukcc", 16)
    params:
        outdir="{sample}/reports/refinement/euk/eukcc",
        db=EUKCC_DB,
        extra=EUKCC_EXTRA,
        eukcc_enabled=EUKCC_ENABLED
    shell:
        r"""
        python - <<'PY'
import sys, shutil, pathlib
eukcc_enabled = "{params.eukcc_enabled}"

if eukcc_enabled != "True":
    sys.stderr.write("ERROR: refinement.euk.eukcc is false, but eukcc_euk_candidates was requested.\n")
    sys.exit(1)

indir = pathlib.Path(r"{input.euk_dir}")
outdir = pathlib.Path(r"{params.outdir}")
done = pathlib.Path(r"{output.done}")
summary = pathlib.Path(r"{output.summary}")

if outdir.exists():
    shutil.rmtree(outdir)
outdir.mkdir(parents=True, exist_ok=True)

fasta_exts = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
files = []
if indir.exists():
    for p in sorted(indir.iterdir()):
        name = p.name.lower()
        if p.is_file() and any(name.endswith(ext) for ext in fasta_exts):
            files.append(p)

if len(files) == 0:
    with open(summary, "w") as out:
        out.write("bin\tcompleteness\tcontamination\twarning\n")
    done.touch()
    sys.exit(0)
PY

        if [ ! -d "{input.euk_dir}" ] || [ -z "$(find {input.euk_dir} -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' -o -name '*.fa.gz' -o -name '*.fna.gz' -o -name '*.fasta.gz' \) -print -quit)" ]; then
            touch {output.done}
            exit 0
        fi

        if [ -n "{params.db}" ]; then
            export EUKCC2_DB="{params.db}"
        fi

        eukcc folder \
          --out {params.outdir} \
          --threads {threads} \
          {params.extra} \
          {input.euk_dir}

        touch {output.done}
        """

# -------------------------------------------------------------------------
# Select best representative per dRep cluster using EukCC
# -------------------------------------------------------------------------
rule select_best_euk_bins:
    input:
        euk_dir="{sample}/reports/refinement/euk/euk_bins",
        kept_manifest="{sample}/reports/refinement/euk/euk_bins/kept_bins.tsv",
        cluster_info="{sample}/reports/refinement/euk/drep/data_tables/Cdb.csv",
        drep_done="{sample}/reports/refinement/euk/drep/.done",
        eukcc="{sample}/reports/refinement/euk/eukcc/eukcc.csv",
        eukcc_done="{sample}/reports/refinement/euk/eukcc/.done"
    output:
        selected_dir=directory("{sample}/reports/refinement/euk/final_bins/selected_bins"),
        selected_manifest="{sample}/reports/refinement/euk/final_bins/selected_bins.tsv",
        done="{sample}/reports/refinement/euk/final_bins/.done"
    shell:
        r"""
        rm -rf {output.selected_dir}
        mkdir -p {output.selected_dir}

        python {SCRIPTS_DIR}/select_best_euk_bins.py \
          --input-dir {input.euk_dir} \
          --clusters {input.cluster_info} \
          --eukcc {input.eukcc} \
          --output-dir {output.selected_dir} \
          --manifest {output.selected_manifest}

        touch {output.done}
        """

# -------------------------------------------------------------------------
# BAT/CAT_pack on final selected bins
# -------------------------------------------------------------------------
rule bat_classify:
    container:
        CONTAINERS["euk"]
    input:
        final_done="{sample}/reports/refinement/euk/final_bins/.done",
        bins_dir="{sample}/reports/refinement/euk/final_bins/selected_bins"
    output:
        classification="{sample}/reports/refinement/euk/bat/bin2classification.txt",
        done="{sample}/reports/refinement/euk/bat/.done"
    threads:
        P["threads"].get("bat", 24)
    params:
        outdir="{sample}/reports/refinement/euk/bat",
        db=BAT_DB,
        taxonomy=BAT_TAXONOMY,
        extra=BAT_EXTRA,
        suffix=BAT_BIN_SUFFIX,
        bat_enabled=BAT_ENABLED,
        db_set=bool(BAT_DB),
        tax_set=bool(BAT_TAXONOMY),
        bat_exe=BAT_EXECUTABLE
    shell:
        r"""
        python - <<'PY'
import sys, shutil, pathlib
bat_enabled = "{params.bat_enabled}"
db_set = "{params.db_set}"
tax_set = "{params.tax_set}"

if bat_enabled != "True":
    sys.stderr.write("ERROR: refinement.euk.bat is false, but bat_classify was requested.\n")
    sys.exit(1)

outdir = pathlib.Path(r"{params.outdir}")
if outdir.exists():
    shutil.rmtree(outdir)
outdir.mkdir(parents=True, exist_ok=True)

if db_set != "True":
    sys.stderr.write("ERROR: refinement.euk.bat_params.db_path is empty.\n")
    sys.exit(1)

if tax_set != "True":
    sys.stderr.write("ERROR: refinement.euk.bat_params.taxonomy_path is empty.\n")
    sys.exit(1)
PY

        if [ ! -d "{input.bins_dir}" ] || [ -z "$(find "{input.bins_dir}" -maxdepth 1 -type f -name "*{params.suffix}" -print -quit)" ]; then
            printf "Bin\tclassification\treason\n" > "{output.classification}"
            touch "{output.done}"
            exit 0
        fi

        {params.bat_exe} bins \
          -b "{input.bins_dir}" \
          -d "{params.db}" \
          -t "{params.taxonomy}" \
          -o "{params.outdir}/bat" \
          -s "{params.suffix}" \
          -n {threads} \
          {params.extra}

        cp "{params.outdir}/bat.bin2classification.txt" "{output.classification}"
        touch "{output.done}"
        """
