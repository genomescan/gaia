# -------------------------------------------------------------------------
# Mapping reads back to metaFlye assembly
# -------------------------------------------------------------------------
rule map_reads_to_assembly:
    container:
        CONTAINERS["metabat2"]
    input:
        assembly=assembly_path("MetaFlye", "{sample}", "assembly.fasta"),
        reads=downstream_reads
    output:
        bam=alignment_path("{sample}", "{sample}.vs_assembly.sorted.bam"),
        bai=alignment_path("{sample}", "{sample}.vs_assembly.sorted.bam.bai")
    threads:
        P["threads"].get("minimap2", 24)
    params:
        preset=MINIMAP2_PRESET,
        extra=MINIMAP2_EXTRA,
        mapping_enabled=MAPPING_ENABLED,
        samtools_threads=max(1, P["threads"].get("samtools", 8))
    shell:
        r"""
        python - <<'PY'
import sys
mapping_enabled = "{params.mapping_enabled}"
if mapping_enabled != "True":
    sys.stderr.write("ERROR: mapping.enabled is false, but map_reads_to_assembly was requested.\n")
    sys.exit(1)
PY
        mkdir -p $(dirname {output.bam})
        minimap2 -ax {params.preset} -t {threads} {params.extra} {input.assembly} {input.reads} \
          | samtools sort -@ {params.samtools_threads} -o {output.bam} -
        samtools index {output.bam}
        """

# -------------------------------------------------------------------------
# Contig depth for downstream binning
# -------------------------------------------------------------------------
rule contig_depth:
    container:
        CONTAINERS["metabat2"]
    input:
        bam=alignment_path("{sample}", "{sample}.vs_assembly.sorted.bam"),
        bai=alignment_path("{sample}", "{sample}.vs_assembly.sorted.bam.bai")
    output:
        depth=alignment_path("{sample}", "{sample}.depth.txt")
    threads:
        P["threads"].get("depth", 8)
    params:
        extra=DEPTH_EXTRA,
        mapping_enabled=MAPPING_ENABLED
    shell:
        r"""
        python - <<'PY'
import sys
mapping_enabled = "{params.mapping_enabled}"
if mapping_enabled != "True":
    sys.stderr.write("ERROR: mapping.enabled is false, but contig_depth was requested.\n")
    sys.exit(1)
PY
        mkdir -p $(dirname {output.depth})
        jgi_summarize_bam_contig_depths \
          --outputDepth {output.depth} \
          {params.extra} \
          {input.bam}
        """
