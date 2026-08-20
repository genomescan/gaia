# -------------------------------------------------------------------------
# Taxonomy: Kraken2
# -------------------------------------------------------------------------
rule kraken2_classify:
    container:
        CONTAINERS["taxonomy"]
    input:
        reads=downstream_reads
    output:
        out=taxonomy_path("Kraken2", "{sample}", "{sample}.kraken2.output.txt"),
        report=taxonomy_path("Kraken2", "{sample}", "{sample}.kraken2.report.txt")
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

        mkdir -p "$(dirname {output.out})"
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
        out=taxonomy_path("Centrifuger", "{sample}", "{sample}.centrifuger.classification.tsv")
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

        mkdir -p "$(dirname {output.out})"
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
        classif=taxonomy_path("Centrifuger", "{sample}", "{sample}.centrifuger.classification.tsv")
    output:
        report=taxonomy_path("Centrifuger", "{sample}", "{sample}.centrifuger.report.tsv")
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
