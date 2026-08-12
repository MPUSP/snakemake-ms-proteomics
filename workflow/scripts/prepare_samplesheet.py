#!/usr/bin/python3

# PREPARE SAMPLESHEET
# -----------------------------------------------------------------------------
#
# This script imports and parses the sample sheeet.
# These tasks include mainly to check that the user-supplied format
# and options are in agreement with the pipeline requirements

from pathlib import Path
import pandas as pd
import re

input_path = snakemake.input["path"]
output_path = snakemake.output["path"]
output_log = snakemake.log["path"]
log = []


def fail(message: str) -> None:
    raise ValueError(message)


def replace_symbols(value: str) -> str:
    return re.sub("[;:,. -]", "_", value)


# import samplesheet
if not Path(input_path).exists():
    fail(
        "Sample sheet was not found under the given path.\n"
        "Please provide a valid path to the file."
    )

df: pd.DataFrame = pd.read_csv(input_path, delimiter="\t")
if len(df.columns) == 1:
    # try import as csv file when table has only 1 col
    df = pd.read_csv(input_path, delimiter=",")
    if len(df.columns) != 6:
        fail("Table has not the correct delimiter; use tabs or commas")
    log += ["Using commas as the default delimiter"]

if len(df.columns) != 6:
    fail(f"Sample sheet has {len(df.columns)} columns, but 6 are required.")

log += [f"Imported TSV file: {input_path}"]
# check that all columns have correct order
expected_columns = [
    "sample",
    "raw_file",
    "condition",
    "replicate",
    "method",
    "comparison",
]
if list(df.columns) != expected_columns:
    fail(f"Columns are not in the correct order. Expected: {expected_columns}")

# replace all special characters by underscores in condtion and contrast names
df[["condition", "comparison"]] = df[["condition", "comparison"]].map(
    lambda x: replace_symbols(x)
)

# checking properties
log += [f"Import from {input_path} successfull. Checking properties..."]
df["raw_file"] = [str(Path(i).resolve()) for i in df["raw_file"]]
log += ["Converted file paths to absolute paths."]
cond = df["condition"].unique()
log += ["Found {0} conditions: {1}".format(len(cond), ", ".join(cond))]
repl = [len(i) for i in df.groupby("condition")["raw_file"].unique()]
log += [
    "Found minimally {0} and maximally {1} replicates per condition".format(
        min(repl), max(repl)
    )
]
types = ["DDA", "DIA", "GPF-DIA", "DIA-Quant", "DIA-Lib"]
sample_type: dict = df.groupby("method").size().to_dict()
sample_type_str: str = ", ".join(sample_type.keys())
log += [f"Found the following sample types: {sample_type_str}"]
if not all([i in types for i in sample_type.keys()]):
    fail(f"Not all sample types are one of {types}")

# check if all comparisons are possible
if not all([i in df["condition"].to_list() for i in df["comparison"]]):
    fail("Some comparisons have controls that don't appear as sample")
log += ["All comparisons are possible, as all controls appear as sample"]

# export modified sample sheet
df.iloc[:, 1:].to_csv(output_path, header=False, index=False, sep="\t")
log += [f"Wrote samplesheet to: {output_path}."]
log += [f"Module finished successfully"]
log = ["SAMPLESHEET: " + i for i in log]
with open(output_log, "w") as log_file:
    log_file.write("\n".join(log))
