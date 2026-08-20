# -------------------------------------------------------------------------
# Raw QC
# -------------------------------------------------------------------------
rule nanoplot_raw:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=raw_fastq
    output:
        html=qc_path("Raw", "{sample}", "{sample}_NanoPlot-report.html"),
        stats=qc_path("Raw", "{sample}", "{sample}_NanoStats.txt")
    threads:
        P["threads"]["nanoplot"]
    shell:
        r"""
        mkdir -p "$(dirname {output.html})"
        NanoPlot \
            --fastq {input.fastq} \
            -o "$(dirname {output.html})" \
            --prefix {wildcards.sample}_
        """


rule nanoqc_raw:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=raw_fastq
    output:
        html=qc_path("Raw", "{sample}", "{sample}_nanoQC.html")
    threads:
        P["threads"]["nanoqc"]
    shell:
        r"""
        mkdir -p "$(dirname {output.html})"
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
        html=qc_path("Filtered", "{sample}", "{sample}_NanoPlot-report.html"),
        stats=qc_path("Filtered", "{sample}", "{sample}_NanoStats.txt")
    threads:
        P["threads"]["nanoplot"]
    shell:
        r"""
        mkdir -p "$(dirname {output.html})"
        NanoPlot \
            --fastq {input.fastq} \
            -o "$(dirname {output.html})" \
            --prefix {wildcards.sample}_
        """


rule nanoqc_filtered:
    container:
        CONTAINERS["preprocessing_qc"]
    input:
        fastq=downstream_reads
    output:
        html=qc_path("Filtered", "{sample}", "{sample}_nanoQC.html")
    threads:
        P["threads"]["nanoqc"]
    shell:
        r"""
        mkdir -p "$(dirname {output.html})"
        tmpdir=$(mktemp -d)

        nanoQC \
            -o "$tmpdir" \
            {input.fastq}

        mv "$tmpdir"/nanoQC.html {output.html}
        rm -rf "$tmpdir"
        """
