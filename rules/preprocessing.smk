# -------------------------------------------------------------------------
# Host removal (minimap2 wrapper)
# -------------------------------------------------------------------------
rule host_removal:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        host_ref=HOST_REF,
        reads=raw_fastq
    output:
        bam="outputs/{sample}/reports/preprocessing/{sample}.host.bam",
        flagstat="outputs/{sample}/reports/preprocessing/{sample}.host.flagstat.txt",
        nohost="outputs/{sample}/reports/preprocessing/{sample}.nohost.fastq.gz"
    threads:
        P["threads"]["host_removal"]
    params:
        out_prefix="outputs/{sample}/reports/preprocessing/{sample}"
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/preprocessing
        scripts/host_removal_mm2.sh \
            -r {input.host_ref} \
            -i {input.reads} \
            -o {params.out_prefix} \
            -t {threads}
        """


# -------------------------------------------------------------------------
# Chopper: length/quality filtering
# -------------------------------------------------------------------------
rule chopper:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        reads="outputs/{sample}/reports/preprocessing/{sample}.nohost.fastq.gz"
    output:
        "outputs/{sample}/reports/preprocessing/{sample}.chopper.fastq.gz"
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
# Filtlong: additional filtering/subsampling
# -------------------------------------------------------------------------
rule filtlong:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        "outputs/{sample}/reports/preprocessing/{sample}.chopper.fastq.gz"
    output:
        "outputs/{sample}/reports/preprocessing/{sample}-preprocessed.fastq.gz"
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
            {input} \
        | gzip -c > {output}
        """
