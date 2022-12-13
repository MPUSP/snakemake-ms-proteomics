#!/usr/bin/python3

# PREPARE DATABASE
# -----------------------------------------------------------------------------
#
# This script attempts to download genome/proteome fasta files
# from NCBI using the NCBI datasets API.
# I a file is provided, nothing happens except that the script
# checks if it is a valid protein FASTA file

from os import path
from io import StringIO
import pandas as pd
import subprocess


input_term = snakemake.params["term"]
output_path = snakemake.output["path"]


if not path.exists(input_term):
    ncbi_result = subprocess.getoutput(
        f"datasets summary genome accession {input_term} --as-json-lines | "
        + "dataformat tsv genome --fields accession,annotinfo-release-date,organism-name"
    )
    if ncbi_result.startswith("Error"):
        print("DATABASE: " + ncbi_result)
        print(
            "DATABASE: the supplied refseq/genbank ID was not valid. Try for example 'GCF_000009045.1'."
        )
        raise ValueError()
    else:
        ncbi_genome = [i.split("\t") for i in ncbi_result.split("\n")]
        ncbi_genome = dict(zip(ncbi_genome[0], ncbi_genome[1]))
        print("DATABASE: Found the following genome(s):\n")
        for k in ncbi_genome.keys():
            print("{0}: {1}".format(k, ncbi_genome.get(k)))
        refseq_id = ncbi_genome.get("Assembly Accession")
        ncbi_command = (
            f"datasets download genome accession {refseq_id}"
            + f" --filename {output_path}/database.zip --include protein; "
            + f"cd {output_path}; unzip database.zip; rm database.zip; "
            + f"cp ncbi_dataset/data/{refseq_id}/protein.faa database.fasta"
        )
        str_out = subprocess.getoutput(ncbi_command)
else:
    # import fasta file
    with open(input_term, "r") as fasta_file:
        fasta = fasta_file.read()

    # check fasta file
    n_items = fasta.count(">")
    print(f"DATABASE: Supplied fasta file '{input_term}' was found.")
    print(f"DATABASE: Supplied fasta file contains {n_items} protein entries.")
    if fasta.count(">rev_"):
        print(
            "DATABASE: Supplied fasta file seems to contain decoy "
            + "proteins, prefix: 'rev_'. Adding decoys is omitted."
        )

    # export fasta file
    with open(path.join(output_path, "database.fasta"), "w") as fasta_out:
        fasta_out.write(fasta)

# print error/log messages
print(f"DATABASE: Module finished successfully")
