# -------------------------------------------------------------------------
# Raw QC: NanoPlot + NanoQC
# -------------------------------------------------------------------------
rule nanoplot_raw:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=raw_fastq
    output:
        html="outputs/{sample}/reports/QC/{sample}_RAW_NanoPlot-report.html"
    threads:
        P["threads"]["nanoplot"]
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/QC
        NanoPlot \
            --fastq {input.fastq} \
            -o outputs/{wildcards.sample}/reports/QC \
            --prefix {wildcards.sample}_RAW

        mv outputs/{wildcards.sample}/reports/QC/{wildcards.sample}_RAWNanoPlot-report.html {output.html}
        """


rule nanoqc_raw:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=raw_fastq
    output:
        html="outputs/{sample}/reports/QC/{sample}_RAW_nanoQC.html"
    threads:
        P["threads"]["nanoqc"]
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/QC
        tmpdir=$(mktemp -d)

        nanoQC \
            -o "$tmpdir" \
            {input.fastq}

        mv "$tmpdir"/nanoQC.html {output.html}
        rm -rf "$tmpdir"
        """


# -------------------------------------------------------------------------
# Filtered (downstream) QC: NanoPlot + NanoQC
# -------------------------------------------------------------------------
rule nanoplot_filtered:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=downstream_reads
    output:
        html="outputs/{sample}/reports/QC/{sample}_NanoPlot-report.html"
    threads:
        P["threads"]["nanoplot"]
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/QC
        NanoPlot \
            --fastq {input.fastq} \
            -o outputs/{wildcards.sample}/reports/QC \
            --prefix {wildcards.sample}

        mv outputs/{wildcards.sample}/reports/QC/{wildcards.sample}NanoPlot-report.html {output.html}
        """


rule nanoqc_filtered:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=downstream_reads
    output:
        html="outputs/{sample}/reports/QC/{sample}_nanoQC.html"
    threads:
        P["threads"]["nanoqc"]
    shell:
        r"""
        mkdir -p outputs/{wildcards.sample}/reports/QC
        tmpdir=$(mktemp -d)

        nanoQC \
            -o "$tmpdir" \
            {input.fastq}

        mv "$tmpdir"/nanoQC.html {output.html}
        rm -rf "$tmpdir"
        """
