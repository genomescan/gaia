# -------------------------------------------------------------------------
# BAM → FASTQ conversion (only triggered when input is a .bam file)
# -------------------------------------------------------------------------
rule convert_bam_to_fastq:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        bam=lambda wc: _search_input_file(wc.sample)
    output:
        fastq="{sample}/reports/preprocessing/{sample}.input.fastq.gz"
    threads:
        P["threads"]["samtools"]
    shell:
        r"""
        mkdir -p {wildcards.sample}/reports/preprocessing
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
            bam="{sample}/reports/preprocessing/{sample}.host.bam",
            flagstat="{sample}/reports/preprocessing/{sample}.host.flagstat.txt",
            nohost="{sample}/reports/preprocessing/{sample}.nohost.fastq.gz"
        threads:
            P["threads"]["host_removal"]
        params:
            out_prefix="{sample}/reports/preprocessing/{sample}"
        shell:
            r"""
            mkdir -p {wildcards.sample}/reports/preprocessing
            {SCRIPTS_DIR}/host_removal_mm2.sh \
                -r {input.host_ref} \
                -i {input.reads} \
                -o {params.out_prefix} \
                -t {threads}
            """

    rule parse_host_removal_stats:
        input:
            flagstat="{sample}/reports/preprocessing/{sample}.host.flagstat.txt"
        output:
            json="{sample}/reports/preprocessing/{sample}.host_removal_stats.json"
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
        return f"{wc.sample}/reports/preprocessing/{wc.sample}.nohost.fastq.gz"
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
        "{sample}/reports/preprocessing/{sample}.chopper.fastq.gz"
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
        "{sample}/reports/preprocessing/{sample}.filtlong.fastq.gz"
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


# -------------------------------------------------------------------------
# Final preprocessed reads: symlink/copy from selected filtering method
# -------------------------------------------------------------------------
rule finalize_preprocessed:
    input:
        reads=lambda wc: (
            f"{wc.sample}/reports/preprocessing/{wc.sample}.chopper.fastq.gz"
            if FILTERING_METHOD == "chopper"
            else f"{wc.sample}/reports/preprocessing/{wc.sample}.filtlong.fastq.gz"
        )
    output:
        "{sample}/reports/preprocessing/{sample}-preprocessed.fastq.gz"
    shell:
        r"""
        cp {input.reads} {output}
        """
