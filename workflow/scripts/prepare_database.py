#!/usr/bin/python3

from pathlib import Path
from Bio import SeqIO

input_fasta = snakemake.input["fasta"]
output_fasta = snakemake.output["fasta"]
output_log = snakemake.log["path"]
log = []
error = []


# read the provided FASTA file
with open(input_fasta) as handle:
    records = [r for r in SeqIO.parse(handle, "fasta")]
log += [f"Supplied fasta file '{input_fasta}' was found"]

# check basic stats
n_items = len(records)
if n_items:
    log += [f"Supplied fasta file contains {n_items} protein entries"]
    decoy_prefix = ["XXX_", "rev_", "Rev_", "REV_"]
    for prefix in decoy_prefix:
        if any(r.id.startswith(prefix) for r in records):
            log += [
                "Supplied fasta file seems to contain decoy "
                + f"proteins with prefix: '>{prefix}'. Adding decoys is omitted"
            ]
            if prefix != "rev_":
                for r in records:
                    if r.id.startswith(prefix):
                        r.id = r.id.replace(prefix, "rev_")
                        log += [
                            f"Replaced decoy prefix '{prefix}' with standard prefix '>rev_' for record '{r.id}'"
                        ]
    if all(not r.id.startswith(p) for p in decoy_prefix for r in records):
        log += [
            f"File does not contain any of the decoy prefixes '{', '.join(decoy_prefix)}'",
            "Decoys will be added by 'decoypyrat'",
        ]
else:
    error += ["The supplied fasta file contains no valid entries"]

# export fasta file
with open(output_fasta, "w") as fasta_out:
    SeqIO.write(records, fasta_out, "fasta")

# print error/log messages
if error:
    print("\n".join(error))
    raise ValueError(
        "Location or format of the supplied database entry was not correct, quitting"
    )
else:
    log += [f"Module finished successfully"]
    log = ["DATABASE: " + i for i in log]
    with open(output_log, "w") as log_file:
        log_file.write("\n".join(log))
