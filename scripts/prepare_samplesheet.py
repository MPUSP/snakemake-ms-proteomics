#!/usr/bin/python3

# PREPARE SAMPLESHEET
# -----------------------------------------------------------------------------
#
# This script imports and parses the sample sheeet.
# These tasks include mainly to check that the user-supplied format
# and options are in agreement with the pipeline requirements

from os import path
import pandas as pd
import re


input_path = snakemake.input["path"]
output_path = snakemake.output["path"]
output_log = snakemake.log["path"]
log = []
error = []


if not path.exists(input_path):
    error += [
        """
        Sample sheet was not found under the given path.
        Please provide a valid path to the file.
        """
    ]
else:
    # import samplesheet
    df_sheet = pd.read_csv(input_path, header=None, delimiter="\t")
    if len(df_sheet.columns) == 1:
        # try import as csv file when table has only 1 col
        df_sheet = pd.read_csv(input_path, header=None, delimiter=",")
        if len(df_sheet.columns) == 1:
            error += ["Table has not the correct delimiter; Use tabs or commas"]
        else:
            log += ["Using commas as the default delimiter"]
    # replace all special characters by underscores
    def replace_symbols(x):
        return(re.sub("[;:,. -]", "_", x))
    df_sheet[[1]] = df_sheet[[1]].applymap(lambda x: replace_symbols(x))
    # checking properties
    if len(df_sheet.columns) >= 4:
        log += [f"Import from {input_path} successfull. Checking properties..."]
        df_sheet[0] = [path.abspath(i) for i in df_sheet[0]]
        log += ["Converted file paths to absolute paths."]
        cond = df_sheet[1].unique()
        log += ["Found {0} conditions: {1}".format(len(cond), ", ".join(cond))]
        repl = [len(i) for i in df_sheet.groupby(1)[2].unique()]
        log += [
            "Found minimally {0} and maximally {1} replicates per condition".format(
                min(repl), max(repl)
            )
        ]
        types = ["DDA", "DIA", "GPF-DIA", "DIA-Quant", "DIA-Lib"]
        sample_type = df_sheet.groupby(3).size().to_dict()
        sample_type_str = ", ".join(sample_type.keys())
        log += [f"Found the following sample types: {sample_type_str}"]
        if not all([i in types for i in sample_type.keys()]):
            error += [f"Not all sample types are one of {types}"]

        if len(df_sheet.columns) >= 5:
            df_sheet[[4]] = df_sheet[[4]].applymap(lambda x: replace_symbols(x))
            if all([i in df_sheet[1].to_list() for i in df_sheet[4]]):
                log += ["Column with constrasts detected"]
            else:
                error += ["Some constrasts have controls that don't appear as sample"]

        # export modified sample sheet
        df_sheet.to_csv(output_path, header=None, index=None, sep="\t")
        log += [f"Wrote samplesheet to: {output_path}."]

# print error/log messages
if error:
    print("\n".join(error))
    raise ValueError(
        "Location or formatting of the sample sheet was not correct, quitting."
    )
else:
    log += [f"Module finished successfully"]
    log = ["SAMPLESHEET: " + i for i in log]
    with open(output_log, "w") as log_file:
        log_file.write("\n".join(log))
