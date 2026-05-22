# -------------------------------------------------------------------------
# Taxonomy: Kraken2
# -------------------------------------------------------------------------
rule kraken2_classify:
    container:
        CONTAINERS["taxonomy"]
    input:
        reads=downstream_reads
    output:
        out="outputs/{sample}/reports/taxonomy/kraken2/{sample}.kraken2.output.txt",
        report="outputs/{sample}/reports/taxonomy/kraken2/{sample}.kraken2.report.txt"
    threads:
        P["threads"].get("kraken2", 16)
    params:
        db=KRAKEN2_DB,
        extra=KRAKEN2_EXTRA,
        enabled=("kraken2" in TAX_TOOLS)
    shell:
        r"""
        if [ "{params.enabled}" = "True" ] && [ -z "{params.db}" ]; then
            echo "ERROR: taxonomy.kraken2_db is empty but kraken2 is enabled." >&2
            exit 1
        fi

        mkdir -p outputs/{wildcards.sample}/reports/taxonomy/kraken2
        kraken2 \
          --db {params.db} \
          --threads {threads} \
          --report {output.report} \
          {params.extra} \
          {input.reads} \
          > {output.out}
        """

# -------------------------------------------------------------------------
# Taxonomy: Centrifuger
# -------------------------------------------------------------------------
rule centrifuger_classify:
    container:
        CONTAINERS["taxonomy"]
    input:
        reads=downstream_reads,
	kraken_done=centrifuger_serial_dep
    output:
        out="outputs/{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.classification.tsv"
    threads:
        P["threads"].get("centrifuger", 16)
    params:
        db=CENTRIFUGER_DB,
        extra=CENTRIFUGER_EXTRA,
        enabled=("centrifuger" in TAX_TOOLS)
    shell:
        r"""
        if [ "{params.enabled}" = "True" ] && [ -z "{params.db}" ]; then
            echo "ERROR: taxonomy.centrifuger_db is empty but centrifuger is enabled." >&2
            exit 1
        fi

        mkdir -p outputs/{wildcards.sample}/reports/taxonomy/centrifuger
        centrifuger \
          -x {params.db} \
          -u {input.reads} \
          -t {threads} \
          {params.extra} \
          > {output.out}
        """

rule centrifuger_quant:
    container:
        CONTAINERS["taxonomy"]
    input:
        classif="outputs/{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.classification.tsv"
    output:
        report="outputs/{sample}/reports/taxonomy/centrifuger/{sample}.centrifuger.report.tsv"
    threads:
        P["threads"].get("centrifuger_quant", 4)
    params:
        db=CENTRIFUGER_DB,
        extra=CENTRIFUGER_QUANT_EXTRA,
        enabled=("centrifuger" in TAX_TOOLS)
    shell:
        r"""
        if [ "{params.enabled}" = "True" ] && [ -z "{params.db}" ]; then
            echo "ERROR: taxonomy.centrifuger_db is empty but centrifuger is enabled." >&2
            exit 1
        fi

        centrifuger-quant \
          -x {params.db} \
          -c {input.classif} \
          {params.extra} \
          > {output.report}
        """
