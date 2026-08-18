import os

SAMPLES = config["samples"]
SCRIPTS_DIR = config.get("path_scripts", os.path.join(os.path.dirname(workflow.snakefile), "scripts"))
RAW_DIR = config["paths"]["raw_dir"]
HOST_REF = config["paths"]["host_ref"]
LAMBDA_REF = config["paths"].get("lambda_ref", "")
CONTAINERS = config.get("containers", {})
P = config["params"]

# -------------------------------------------------------------------------
# Run mode config
# -------------------------------------------------------------------------
RUN_CFG = config.get("run", {})
RUN_MODE = RUN_CFG.get("mode", "both")
SERIAL_PROFILE_THEN_ASSEMBLY = RUN_CFG.get("serial_profile_then_assembly", True)

VALID_RUN_MODES = {"profiling_only", "assembly_binning_only", "both"}
if RUN_MODE not in VALID_RUN_MODES:
    raise ValueError(
        f"Invalid run.mode: {RUN_MODE}. "
        f"Choose one of: {', '.join(sorted(VALID_RUN_MODES))}"
    )

RUN_PROFILE = RUN_MODE in {"profiling_only", "both"}
RUN_ASSEMBLY = RUN_MODE in {"assembly_binning_only", "both"}

PREPROCESSING_CFG = config.get("preprocessing", {})
PREPROCESSING_ENABLED = PREPROCESSING_CFG.get("enabled", True)

def raw_fastq(wc):
    gz = os.path.join(RAW_DIR, f"{wc.sample}.fastq.gz")
    plain = os.path.join(RAW_DIR, f"{wc.sample}.fastq")
    if os.path.exists(gz):
        return gz
    if os.path.exists(plain):
        return plain
    return gz

def downstream_reads(wc):
    if PREPROCESSING_ENABLED:
        return f"outputs/{wc.sample}/reports/preprocessing/{wc.sample}-preprocessed.fastq.gz"
    return raw_fastq(wc)

def profile_stage_done(wc):
    if RUN_PROFILE and RUN_ASSEMBLY and SERIAL_PROFILE_THEN_ASSEMBLY:
        return f"outputs/{wc.sample}/stages/profile/{wc.sample}.done"
    return []

# -------------------------------------------------------------------------
# Taxonomy config
# -------------------------------------------------------------------------
TAX = config.get("taxonomy", {})
TAX_TOOLS = set(TAX.get("tools", []))

KRAKEN2_DB = TAX.get("kraken2_db", "")
CENTRIFUGER_DB = TAX.get("centrifuger_db", "")

KRAKEN2_EXTRA = TAX.get("kraken2", {}).get("extra", "")
CENTRIFUGER_EXTRA = TAX.get("centrifuger", {}).get("extra", "")
CENTRIFUGER_QUANT_EXTRA = TAX.get("centrifuger", {}).get("quant", {}).get("extra", "")

TAXONOMY_ENABLED = RUN_PROFILE and len(TAX_TOOLS) > 0

def centrifuger_serial_dep(wc):
    if "kraken2" in TAX_TOOLS and "centrifuger" in TAX_TOOLS:
        return f"outputs/{wc.sample}/reports/taxonomy/kraken2/{wc.sample}.kraken2.report.txt"
    return []

# -------------------------------------------------------------------------
# Assembly config
# -------------------------------------------------------------------------
ASM = config.get("assembly", {})
ASSEMBLY_ENABLED = RUN_ASSEMBLY
ASSEMBLY_TOOL = ASM.get("tool", "metaflye")
ASSEMBLY_READ_TYPE = ASM.get("read_type", "--nano-hq")
ASSEMBLY_EXTRA = ASM.get("extra", "")

ASM_QC = config.get("assembly_qc", {})
METAQUAST_ENABLED = ASSEMBLY_ENABLED and ASM_QC.get("metaquast", False)
METAQUAST_REFERENCES = ASM_QC.get("metaquast_references", "")
METAQUAST_EXTRA = ASM_QC.get("extra", "")

# -------------------------------------------------------------------------
# Mapping / depth config
# -------------------------------------------------------------------------
MAP_CFG = config.get("mapping", {})
MAPPING_ENABLED = RUN_ASSEMBLY and MAP_CFG.get("enabled", ASSEMBLY_ENABLED)
MINIMAP2_PRESET = MAP_CFG.get("preset", "map-ont")
MINIMAP2_EXTRA = MAP_CFG.get("extra", "")
DEPTH_EXTRA = MAP_CFG.get("depth_extra", "")

# -------------------------------------------------------------------------
# Binning config
# -------------------------------------------------------------------------
BIN_CFG = config.get("binning", {})
BINNING_ENABLED = RUN_ASSEMBLY and BIN_CFG.get("enabled", ASSEMBLY_ENABLED and MAPPING_ENABLED)
BINNING_TOOLS = set(BIN_CFG.get("tools", [])) if BINNING_ENABLED else set()

METABAT2_CFG = BIN_CFG.get("metabat2", {})
METABAT2_MIN_CONTIG = METABAT2_CFG.get("min_contig", 1500)
METABAT2_EXTRA = METABAT2_CFG.get("extra", "")

SEMIBIN2_CFG = BIN_CFG.get("semibin2", {})
SEMIBIN2_EXTRA = SEMIBIN2_CFG.get("extra", "")

COMEBIN_CFG = BIN_CFG.get("comebin", {})
COMEBIN_VIEWS = COMEBIN_CFG.get("views", 6)
COMEBIN_EXTRA = COMEBIN_CFG.get("extra", "")

# -------------------------------------------------------------------------
# DAS Tool dynamic helpers
# -------------------------------------------------------------------------
VALID_DASTOOL_BINS = ["metabat2", "semibin2", "comebin"]

def enabled_dastool_tools():
    return [tool for tool in VALID_DASTOOL_BINS if tool in BINNING_TOOLS]

def dastool_inputs(wc):
    return [
        f"outputs/{wc.sample}/reports/binning/normalized/{tool}.contig2bin.tsv"
        for tool in enabled_dastool_tools()
    ]

def dastool_labels():
    return ",".join(enabled_dastool_tools())

def dastool_input_string(wc):
    return ",".join(dastool_inputs(wc))

# -------------------------------------------------------------------------
# Prokaryotic refinement config
# -------------------------------------------------------------------------
REF_CFG = config.get("refinement", {})
PROK_CFG = REF_CFG.get("prok", {})

DASTOOL_ENABLED = RUN_ASSEMBLY and PROK_CFG.get("dastool", BINNING_ENABLED)
CHECKM2_ENABLED = RUN_ASSEMBLY and PROK_CFG.get("checkm2", DASTOOL_ENABLED)

DASTOOL_CFG = PROK_CFG.get("dastool_params", {})
DASTOOL_EXTRA = DASTOOL_CFG.get("extra", "")
DASTOOL_SCORE_THRESHOLD = DASTOOL_CFG.get("score_threshold", 0.5)
DASTOOL_SEARCH_ENGINE = DASTOOL_CFG.get("search_engine", "diamond")
DASTOOL_WRITE_BIN_EVALS = DASTOOL_CFG.get("write_bin_evals", True)
DASTOOL_WRITE_UNBINNED = DASTOOL_CFG.get("write_unbinned", False)
DASTOOL_DB_DIRECTORY = DASTOOL_CFG.get("db_directory", "db")

CHECKM2_CFG = PROK_CFG.get("checkm2_params", {})
CHECKM2_DB = CHECKM2_CFG.get("db_path", "")
CHECKM2_EXTENSION = CHECKM2_CFG.get("extension", "fa")
CHECKM2_EXTRA = CHECKM2_CFG.get("extra", "")

GTDBTK_ENABLED = RUN_ASSEMBLY and PROK_CFG.get("gtdbtk", CHECKM2_ENABLED)
GTDBTK_CFG = PROK_CFG.get("gtdbtk_params", {})
GTDBTK_DB = GTDBTK_CFG.get("db_path", "")
GTDBTK_EXTRA = GTDBTK_CFG.get("extra", "")
GTDBTK_MASH_DB = GTDBTK_CFG.get("mash_db", "")

# -------------------------------------------------------------------------
# Eukaryotic refinement config
# -------------------------------------------------------------------------
EUK_CFG = REF_CFG.get("euk", {})

ACR_ENABLED = RUN_ASSEMBLY and EUK_CFG.get("acr", BINNING_ENABLED)

ACR_CFG = EUK_CFG.get("acr_params", {})
ACR_EXTRA = ACR_CFG.get("extra", "")
ACR_PREFIX = ACR_CFG.get("prefix", "refine")
ACR_MIN_SIZE = ACR_CFG.get("min_size", "500k")
ACR_BYPASS = ACR_CFG.get("bypass", "N")
ACR_FROM_JGI = ACR_CFG.get("from_jgi_cov", "Y")
ACR_RUN_GMESEUK = ACR_CFG.get("run_gmesEuk", "Y")
ACR_TARGET = ACR_CFG.get("target", "Both")
ACR_COMP = ACR_CFG.get("comp", 50)
ACR_CONT = ACR_CFG.get("cont", 10)

DREP_ENABLED = RUN_ASSEMBLY and EUK_CFG.get("drep", True)
DREP_CFG = EUK_CFG.get("drep_params", {})
DREP_EXTRA = DREP_CFG.get("extra", "")

EUKCC_ENABLED = RUN_ASSEMBLY and EUK_CFG.get("eukcc", True)
EUKCC_CFG = EUK_CFG.get("eukcc_params", {})
EUKCC_EXTRA = EUKCC_CFG.get("extra", "")
EUKCC_DB = EUKCC_CFG.get("db_path", "")

BAT_ENABLED = RUN_ASSEMBLY and EUK_CFG.get("bat", True)
BAT_CFG = EUK_CFG.get("bat_params", {})
BAT_EXTRA = BAT_CFG.get("extra", "")
BAT_DB = BAT_CFG.get("db_path", "")
BAT_TAXONOMY = BAT_CFG.get("taxonomy_path", "")
BAT_BIN_SUFFIX = BAT_CFG.get("bin_suffix", ".fa")

# -------------------------------------------------------------------------
# Tool executables
# -------------------------------------------------------------------------
FLYE_EXECUTABLE = ASM.get("executable", "flye")
COMEBIN_EXECUTABLE = COMEBIN_CFG.get("executable", "run_comebin.sh")
BAT_EXECUTABLE = BAT_CFG.get("executable", "CAT_pack")
