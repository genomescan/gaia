# -------------------------------------------------------------------------
# Raw QC
# -------------------------------------------------------------------------
rule nanoplot_raw:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=raw_fastq
    output:
        html="outputs/{sample}/reports/QC/raw/{sample}_NanoPlot-report.html"
    threads:
        P["threads"]["nanoplot"]
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/QC/raw
        NanoPlot \
            --fastq {input.fastq} \
            -o outputs/{wildcards.sample}/reports/QC/raw \
            --prefix {wildcards.sample}_
        """


rule nanoqc_raw:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=raw_fastq
    output:
        html="outputs/{sample}/reports/QC/raw/{sample}_nanoQC.html"
    threads:
        P["threads"]["nanoqc"]
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/QC/raw
        tmpdir=$(mktemp -d)

        nanoQC \
            -o "$tmpdir" \
            {input.fastq}

        mv "$tmpdir"/nanoQC.html {output.html}
        rm -rf "$tmpdir"
        """


# -------------------------------------------------------------------------
# Filtered QC
# -------------------------------------------------------------------------
rule nanoplot_filtered:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=downstream_reads
    output:
        html="outputs/{sample}/reports/QC/filtered/{sample}_NanoPlot-report.html"
    threads:
        P["threads"]["nanoplot"]
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/QC/filtered
        NanoPlot \
            --fastq {input.fastq} \
            -o outputs/{wildcards.sample}/reports/QC/filtered \
            --prefix {wildcards.sample}_
        """


rule nanoqc_filtered:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=downstream_reads
    output:
        html="outputs/{sample}/reports/QC/filtered/{sample}_nanoQC.html"
    threads:
        P["threads"]["nanoqc"]
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/QC/filtered
        tmpdir=$(mktemp -d)

        nanoQC \
            -o "$tmpdir" \
            {input.fastq}

        mv "$tmpdir"/nanoQC.html {output.html}
        rm -rf "$tmpdir"
        """
