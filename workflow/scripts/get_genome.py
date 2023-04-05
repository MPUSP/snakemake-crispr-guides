#!/usr/bin/python3

# GET GENOME
# -----------------------------------------------------------------------------
#
# This script attempts to download genome sequence files and
# genome annotation from NCBI using the NCBI datasets API, or an
# alternative database. If a genbank file is provided,
# the script checks if it has a valid format and prepares the data for
# downstream nodules.

from os import path
from io import StringIO
from subprocess import getoutput

input_term = snakemake.params["term"]
input_database = snakemake.params["database"]
output_path = snakemake.output["path"]
output_log = snakemake.log["path"]
log = []
error = []


if not path.exists(input_term):
    ncbi_result = getoutput(
        f"datasets summary genome accession {input_term} --as-json-lines | "
        + "dataformat tsv genome --fields accession,annotinfo-release-date,organism-name"
    )
    if ncbi_result.startswith("Error"):
        error += [ncbi_result]
        error += [
            "The supplied refseq/genbank ID was not valid. Example for correct input: 'GCF_000009045.1'"
        ]
    else:
        ncbi_genome = [
            i.split("\t")
            for i in ncbi_result.split("\n")
            if not i.startswith("New version")
        ]
        ncbi_genome = dict(zip(ncbi_genome[0], ncbi_genome[1]))
        log += ["Found the following genome(s):"]
        for k in ncbi_genome.keys():
            log += ["{0}: {1}".format(k, ncbi_genome.get(k))]
        refseq_id = ncbi_genome.get("Assembly Accession")
        if not refseq_id.startswith("GCF_"):
            error += ["The RefSeq ID '{0}' has no valid format.".format(refseq_id)]
        ncbi_command = (
            f"datasets download genome accession {refseq_id}"
            + f" --filename {output_path}/database.zip --include genome,gff3; "
            + f"cd {output_path}; unzip database.zip; rm database.zip; "
            + f"cp ncbi_dataset/data/{refseq_id}/*.fna genome.fasta; "
            + f"cp ncbi_dataset/data/{refseq_id}/genomic.gff genome.gff"
        )
        str_out = getoutput(ncbi_command)
else:
    # import fasta file
    with open(input_term, "r") as fasta_file:
        fasta = fasta_file.read()

    # check fasta file
    n_items = fasta.count(">")
    if n_items:
        log += [f"Supplied fasta file '{input_term}' was found"]
        log += [f"Supplied fasta file contains {n_items} protein entries"]
        decoy_prefix = [">XXX_", ">rev_", ">Rev_", ">REV_"]
        for prefix in decoy_prefix:
            if fasta.count(prefix):
                log += [
                    "Supplied fasta file seems to contain decoy "
                    + f"proteins with prefix: '{prefix}'. Adding decoys is omitted"
                ]
                if prefix != ">rev_":
                    fasta.replace(prefix, ">rev_")
                    log += [
                        f"Replaced decoy prefix '{prefix}' with standard prefix '>rev_'"
                    ]
        if all([i not in fasta for i in decoy_prefix]):
            log += [
                "File does not contain any of the decoy prefixes '{0}'".format(
                    "', '".join(decoy_prefix)
                ),
                "Decoys will be added by 'decoypyrat'",
            ]
    else:
        error += ["The supplied fasta file contains no valid entries starting with '>'"]

    # export fasta file
    with open(path.join(output_path, "database.fasta"), "w") as fasta_out:
        fasta_out.write(fasta)

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
