# -------------------------------------------------------------------------
# BAM → FASTQ conversion (only triggered when input is a .bam file)
# -------------------------------------------------------------------------
rule convert_bam_to_fastq:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        bam=lambda wc: _search_input_file(wc.sample)
    output:
        fastq=preprocessing_path("{sample}", "{sample}.input.fastq.gz")
    threads:
        P["threads"]["samtools"]
    shell:
        r"""
        mkdir -p "$(dirname {output.fastq})"
        samtools fastq -@ {threads} {input.bam} -T "*" -0 {output.fastq}
        """


# -------------------------------------------------------------------------
# Host removal (minimap2 wrapper) – only runs when HOST_REF is provided
# -------------------------------------------------------------------------
if HOST_REMOVAL_ENABLED:
    rule host_removal:
        container:
            CONTAINERS["preprocessing_qc"]
        input:
            host_ref=HOST_REF,
            reads=raw_fastq_or_converted
        output:
            bam=preprocessing_path("{sample}", "{sample}.host.bam"),
            flagstat=preprocessing_path("{sample}", "{sample}.host.flagstat.txt"),
            nohost=preprocessing_path("{sample}", "{sample}.nohost.fastq.gz")
        threads:
            P["threads"]["host_removal"]
        params:
            out_prefix=preprocessing_path("{sample}", "{sample}")
        shell:
            r"""
            mkdir -p "$(dirname {output.bam})"
            {SCRIPTS_DIR}/host_removal_mm2.sh \
                -r {input.host_ref} \
                -i {input.reads} \
                -o {params.out_prefix} \
                -t {threads}
            """

    rule parse_host_removal_stats:
        input:
            flagstat=preprocessing_path("{sample}", "{sample}.host.flagstat.txt")
        output:
            json=preprocessing_path("{sample}", "{sample}.host_removal_stats.json")
        shell:
            r"""
            python {SCRIPTS_DIR}/parse_host_removal.py \
                --sample {wildcards.sample} \
                --flagstat {input.flagstat} \
                --output {output.json}
            """


def _reads_after_host_removal(wc):
    """Return the reads to use for filtering (post-host-removal or raw)."""
    if HOST_REMOVAL_ENABLED:
        return preprocessing_path(wc.sample, f"{wc.sample}.nohost.fastq.gz")
    return raw_fastq_or_converted(wc)


# -------------------------------------------------------------------------
# Chopper: length/quality filtering (default method)
# -------------------------------------------------------------------------
rule chopper_filter:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        reads=_reads_after_host_removal
    output:
        preprocessing_path("{sample}", "{sample}.chopper.fastq.gz")
    threads:
        P["threads"]["chopper"]
    params:
        min_length=P["chopper"]["min_length"],
        q_threshold=P["chopper"]["quality_threshold"],
    shell:
        r"""
        gunzip -c {input.reads} \
          | chopper \
              --threads {threads} \
              --minlength {params.min_length} \
              --quality {params.q_threshold} \
          | gzip -c > {output}
        """


# -------------------------------------------------------------------------
# Filtlong: alternative filtering method
# -------------------------------------------------------------------------
rule filtlong_filter:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        reads=_reads_after_host_removal
    output:
        preprocessing_path("{sample}", "{sample}.filtlong.fastq.gz")
    threads:
        P["threads"]["filtlong"]
    params:
        min_length=P["filtlong"]["min_length"],
        keep_percent=P["filtlong"]["keep_percent"]
    shell:
        r"""
        filtlong \
            --min_length {params.min_length} \
            --keep_percent {params.keep_percent} \
            {input.reads} \
        | gzip -c > {output}
        """

