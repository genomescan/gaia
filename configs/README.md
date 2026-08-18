## Configuration files

This directory contains configuration files that define the pipeline workflow variant to run.
Select one with `-c / --config` when calling the `run` wrapper.

Each config file is plain JSON (no extension) and must contain at minimum:

| Key | Description |
|-----|-------------|
| `snakefile` | Path to the Snakemake entry point, relative to the pipeline root (required) |
| `description` | Human-readable description shown in `--help` |
| `run_mode` | Workflow mode: `both`, `profiling_only`, or `assembly_binning_only` |
| `serial_profile_then_assembly` | Whether to run profiling before assembly (`true`/`false`) |
| `taxonomy_tools` | List of taxonomy tools to run, e.g. `["kraken2", "centrifuger"]` |
| `kraken2_db` | Path to Kraken2 database directory |
| `centrifuger_db` | Path to Centrifuger database prefix |
| `binning_tools` | List of binning tools, e.g. `["metabat2", "semibin2", "comebin"]` |

To create a new config, copy an existing one and modify accordingly.
