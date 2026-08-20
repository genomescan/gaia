# -------------------------------------------------------------------------
# Assembly: metaFlye
# -------------------------------------------------------------------------
rule metaflye_assemble:
    container:
        CONTAINERS["assembly"]
    input:
        reads=downstream_reads,
        profile_done=profile_stage_done
    output:
        assembly=assembly_path("MetaFlye", "{sample}", "assembly.fasta"),
        assembly_info=assembly_path("MetaFlye", "{sample}", "assembly_info.txt"),
        graph=assembly_path("MetaFlye", "{sample}", "assembly_graph.gfa")
    threads:
        P["threads"].get("metaflye", 32)
    params:
        outdir=assembly_path("MetaFlye", "{sample}"),
        read_type=ASSEMBLY_READ_TYPE,
        extra=ASSEMBLY_EXTRA,
        enabled=ASSEMBLY_ENABLED,
        flye_exe=FLYE_EXECUTABLE
    shell:
        r"""
        python - <<'PY'
import sys

enabled = "{params.enabled}"

if enabled != "True":
    sys.stderr.write(
        "ERROR: assembly.enabled is false, but rule metaflye_assemble was requested.\n"
    )
    sys.exit(1)
PY

        mkdir -p {params.outdir}

        {params.flye_exe} \
          --meta \
          {params.read_type} \
          {input.reads} \
          --threads {threads} \
          --out-dir {params.outdir} \
          {params.extra}
        """


# -------------------------------------------------------------------------
# Assembly QC: optional MetaQUAST
# -------------------------------------------------------------------------
rule metaquast_assembly:
    container:
        CONTAINERS["assembly"]
    input:
        assembly=assembly_path("MetaFlye", "{sample}", "assembly.fasta")
    output:
        report=assembly_path("MetaQUAST", "{sample}", "report.html"),
        tsv=assembly_path("MetaQUAST", "{sample}", "report.tsv")
    threads:
        P["threads"].get("metaquast", 16)
    params:
        outdir=assembly_path("MetaQUAST", "{sample}"),
        refs=METAQUAST_REFERENCES,
        extra=METAQUAST_EXTRA,
        enabled=METAQUAST_ENABLED
    shell:
        r"""
        python - <<'PY'
import sys

enabled = "{params.enabled}"
refs = "{params.refs}"

if enabled != "True":
    sys.stderr.write(
        "ERROR: assembly_qc.metaquast is false, but rule metaquast_assembly was requested.\n"
    )
    sys.exit(1)

if not refs:
    sys.stderr.write(
        "ERROR: assembly_qc.metaquast_references is empty, but MetaQUAST was requested.\n"
    )
    sys.exit(1)
PY

        mkdir -p {params.outdir}

        metaquast.py \
          {input.assembly} \
          -r {params.refs} \
          -o {params.outdir} \
          --threads {threads} \
          {params.extra}
        """
