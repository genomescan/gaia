import os

configfile: "config.yaml"

include: "rules/common.smk"
include: "rules/qc.smk"
include: "rules/preprocessing.smk"
include: "rules/taxonomy.smk"
include: "rules/assembly.smk"
include: "rules/mapping.smk"
include: "rules/binning.smk"
include: "rules/refine_prok.smk"
include: "rules/refine_euk.smk"
include: "rules/reporting.smk"
include: "rules/targets.smk"
