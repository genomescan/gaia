# -------------------------------------------------------------------------
# MetaBAT2 binning
# -------------------------------------------------------------------------
rule metabat2_bin:
    container:
        CONTAINERS["metabat2"]
    input:
        assembly=assembly_path("MetaFlye", "{sample}", "assembly.fasta"),
        depth=alignment_path("{sample}", "{sample}.depth.txt")
    output:
        bins=directory(binning_path("MetaBAT2", "{sample}", "bins")),
        done=binning_path("MetaBAT2", "{sample}", ".done")
    threads:
        P["threads"].get("metabat2", 16)
    params:
        min_contig=METABAT2_MIN_CONTIG,
        extra=METABAT2_EXTRA,
        out_prefix=binning_path("MetaBAT2", "{sample}", "bins", "bin"),
        mapping_enabled=MAPPING_ENABLED,
        binning_enabled=BINNING_ENABLED,
        tool_enabled=("metabat2" in BINNING_TOOLS)
    shell:
        r"""
        python - <<'PY'
import sys
mapping_enabled = "{params.mapping_enabled}"
binning_enabled = "{params.binning_enabled}"
tool_enabled = "{params.tool_enabled}"

if mapping_enabled != "True":
    sys.stderr.write("ERROR: MetaBAT2 binning requires mapping.enabled: true.\n")
    sys.exit(1)
if binning_enabled != "True":
    sys.stderr.write("ERROR: binning.enabled is false, but metabat2_bin was requested.\n")
    sys.exit(1)
if tool_enabled != "True":
    sys.stderr.write("ERROR: 'metabat2' is not listed in binning.tools, but metabat2_bin was requested.\n")
    sys.exit(1)
PY
        rm -rf {output.bins}
        mkdir -p {output.bins}
        metabat2 \
          -i {input.assembly} \
          -a {input.depth} \
          -o {params.out_prefix} \
          -t {threads} \
          {params.extra}
        touch {output.done}
        """

# -------------------------------------------------------------------------
# SemiBin2 binning
# -------------------------------------------------------------------------
rule semibin2_bin:
    container:
        CONTAINERS["semibin2"]
    input:
        assembly=assembly_path("MetaFlye", "{sample}", "assembly.fasta"),
        bam=alignment_path("{sample}", "{sample}.vs_assembly.sorted.bam")
    output:
        bins=directory(binning_path("SemiBin2", "{sample}", "bins")),
        done=binning_path("SemiBin2", "{sample}", ".done")
    threads:
        P["threads"].get("semibin2", 24)
    params:
        extra=SEMIBIN2_EXTRA,
        outdir=binning_path("SemiBin2", "{sample}"),
        mapping_enabled=MAPPING_ENABLED,
        binning_enabled=BINNING_ENABLED,
        tool_enabled=("semibin2" in BINNING_TOOLS)
    shell:
        r"""
        python - <<'PY'
import sys
mapping_enabled = "{params.mapping_enabled}"
binning_enabled = "{params.binning_enabled}"
tool_enabled = "{params.tool_enabled}"

if mapping_enabled != "True":
    sys.stderr.write("ERROR: SemiBin2 binning requires mapping.enabled: true.\n")
    sys.exit(1)
if binning_enabled != "True":
    sys.stderr.write("ERROR: binning.enabled is false, but semibin2_bin was requested.\n")
    sys.exit(1)
if tool_enabled != "True":
    sys.stderr.write("ERROR: 'semibin2' is not listed in binning.tools, but semibin2_bin was requested.\n")
    sys.exit(1)
PY
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        SemiBin2 single_easy_bin \
          -p {threads} \
          -i {input.assembly} \
          -b {input.bam} \
          -o {params.outdir} \
          --sequencing-type long_read \
          {params.extra}

        mkdir -p {output.bins}

        if [ -d "{params.outdir}/output_bins" ]; then
          cp -a {params.outdir}/output_bins/. {output.bins}/
        elif [ -d "{params.outdir}/bins" ]; then
          cp -a {params.outdir}/bins/. {output.bins}/
        else
          echo "ERROR: SemiBin2 finished but no output bin directory was found in {params.outdir}" >&2
          echo "Contents of {params.outdir}:" >&2
          ls -lah {params.outdir} >&2 || true
          exit 1
        fi

        touch {output.done}
        """

# -------------------------------------------------------------------------
# COMEBin binning
# -------------------------------------------------------------------------
rule comebin_bin:
    container:
        CONTAINERS["comebin"]
    input:
        assembly=assembly_path("MetaFlye", "{sample}", "assembly.fasta"),
        bam=alignment_path("{sample}", "{sample}.vs_assembly.sorted.bam")
    output:
        bins=directory(binning_path("COMEBin", "{sample}", "bins")),
        tsv=binning_path("COMEBin", "{sample}", "comebin_res.tsv"),
        done=binning_path("COMEBin", "{sample}", ".done")
    threads:
        P["threads"].get("comebin", 24)
    params:
        extra=COMEBIN_EXTRA,
        views=COMEBIN_VIEWS,
        outdir=binning_path("COMEBin", "{sample}"),
        bamdir=binning_path("COMEBin", "{sample}", "bamfiles"),
        comebin_exe=COMEBIN_EXECUTABLE,
        mapping_enabled=MAPPING_ENABLED,
        binning_enabled=BINNING_ENABLED,
        tool_enabled=("comebin" in BINNING_TOOLS)
    shell:
        r"""
        python - <<'PY'
import sys
mapping_enabled = "{params.mapping_enabled}"
binning_enabled = "{params.binning_enabled}"
tool_enabled = "{params.tool_enabled}"

if mapping_enabled != "True":
    sys.stderr.write("ERROR: COMEBin binning requires mapping.enabled: true.\n")
    sys.exit(1)
if binning_enabled != "True":
    sys.stderr.write("ERROR: binning.enabled is false, but comebin_bin was requested.\n")
    sys.exit(1)
if tool_enabled != "True":
    sys.stderr.write("ERROR: 'comebin' is not listed in binning.tools, but comebin_bin was requested.\n")
    sys.exit(1)
PY
        rm -rf {params.outdir}
        mkdir -p {params.bamdir}
        ln -sf "$(realpath {input.bam})" {params.bamdir}/$(basename {input.bam})

        {params.comebin_exe} \
          -a {input.assembly} \
          -o {params.outdir} \
          -p {params.bamdir} \
          -n {params.views} \
          -t {threads} \
          {params.extra}

        mkdir -p {output.bins}
        if [ -d "{params.outdir}/comebin_res/comebin_res_bins" ]; then
          cp -a {params.outdir}/comebin_res/comebin_res_bins/. {output.bins}/
        fi
        if [ -f "{params.outdir}/comebin_res/comebin_res.tsv" ]; then
          cp {params.outdir}/comebin_res/comebin_res.tsv {output.tsv}
        else
          touch {output.tsv}
        fi
        touch {output.done}
        """

# -------------------------------------------------------------------------
# Collect all raw bins into one pooled directory
# -------------------------------------------------------------------------
rule collect_all_bins:
    input:
        metabat2_done=binning_path("MetaBAT2", "{sample}", ".done") if "metabat2" in BINNING_TOOLS else [],
        semibin2_done=binning_path("SemiBin2", "{sample}", ".done") if "semibin2" in BINNING_TOOLS else [],
        comebin_done=binning_path("COMEBin", "{sample}", ".done") if "comebin" in BINNING_TOOLS else []
    output:
        bins=directory(binning_path("AllBins", "{sample}")),
        manifest=binning_path("AllBins", "{sample}", "manifest.tsv"),
        done=binning_path("AllBins", "{sample}", ".done")
    params:
        metabat2_bins=lambda wc: binning_path("MetaBAT2", wc.sample, "bins") if "metabat2" in BINNING_TOOLS else "",
        semibin2_bins=lambda wc: binning_path("SemiBin2", wc.sample, "bins") if "semibin2" in BINNING_TOOLS else "",
        comebin_bins=lambda wc: binning_path("COMEBin", wc.sample, "bins") if "comebin" in BINNING_TOOLS else "",
        mapping_enabled=MAPPING_ENABLED,
        binning_enabled=BINNING_ENABLED,
        tools_count=len(BINNING_TOOLS)
    shell:
        r"""
        python - <<'PY'
import sys
mapping_enabled = "{params.mapping_enabled}"
binning_enabled = "{params.binning_enabled}"
tools_count = int("{params.tools_count}")

if mapping_enabled != "True":
    sys.stderr.write("ERROR: collect_all_bins requires mapping.enabled: true.\n")
    sys.exit(1)
if binning_enabled != "True":
    sys.stderr.write("ERROR: collect_all_bins requires binning.enabled: true.\n")
    sys.exit(1)
if tools_count == 0:
    sys.stderr.write("ERROR: collect_all_bins requires at least one tool in binning.tools.\n")
    sys.exit(1)
PY
        python {SCRIPTS_DIR}/collect_all_bins.py \
          --sample {wildcards.sample} \
          --outdir {output.bins} \
          --manifest {output.manifest} \
          --metabat2 "{params.metabat2_bins}" \
          --semibin2 "{params.semibin2_bins}" \
          --comebin "{params.comebin_bins}"
        touch {output.done}
        """

# -------------------------------------------------------------------------
# Convert raw binner FASTA outputs to DAS Tool-style contig2bin tables
# -------------------------------------------------------------------------
rule metabat2_contig2bin:
    input:
        bins=binning_path("MetaBAT2", "{sample}", "bins"),
        done=binning_path("MetaBAT2", "{sample}", ".done")
    output:
        tsv=binning_path("Normalized", "{sample}", "metabat2.contig2bin.tsv")
    shell:
        r"""
        python {SCRIPTS_DIR}/bins_to_contig2bin.py \
          --bins {input.bins} \
          --output {output.tsv} \
          --tool metabat2
        """


rule semibin2_contig2bin:
    input:
        bins=binning_path("SemiBin2", "{sample}", "bins"),
        done=binning_path("SemiBin2", "{sample}", ".done")
    output:
        tsv=binning_path("Normalized", "{sample}", "semibin2.contig2bin.tsv")
    shell:
        r"""
        python {SCRIPTS_DIR}/bins_to_contig2bin.py \
          --bins {input.bins} \
          --output {output.tsv} \
          --tool semibin2
        """


rule comebin_contig2bin:
    input:
        bins=binning_path("COMEBin", "{sample}", "bins"),
        done=binning_path("COMEBin", "{sample}", ".done")
    output:
        tsv=binning_path("Normalized", "{sample}", "comebin.contig2bin.tsv")
    shell:
        r"""
        python {SCRIPTS_DIR}/bins_to_contig2bin.py \
          --bins {input.bins} \
          --output {output.tsv} \
          --tool comebin
        """
