#!/usr/bin/python3

# GET GENOME
# -----------------------------------------------------------------------------
#
# This script attempts to download genome sequence (FASTA) and
# genome annotation (GFF / GTF) files from NCBI using the NCBI datasets
# API, or a similar database. Alternatively, a FASTA and GFF file can be
# supplied by the user. Input is roughly checked for validity.

from os import path
from io import StringIO
from subprocess import getoutput

input_database = snakemake.params["database"]
input_assembly = snakemake.params["assembly"]
input_fasta = snakemake.params["fasta"]
input_gff = snakemake.params["gff"]
output_path = snakemake.output["path"]
output_fasta = snakemake.output["fasta"]
output_gff = snakemake.output["gff"]
output_log = snakemake.log["path"]
log = []
error = []


def check_fasta(input_fasta, log=[], error=[]):
    with open(input_fasta, "r") as fasta_file:
        fasta = fasta_file.read()
    n_items = fasta.count(">")
    if n_items:
        log += [f"Supplied fasta file '{input_fasta}' was found"]
        log += [f"Supplied fasta file contains {n_items} items"]
    else:
        error += ["The supplied fasta file contains no valid entries starting with '>'"]
    return fasta, log, error


def check_gff(input_gff, log=[], error=[]):
    with open(input_gff, "r") as gff_file:
        gff = gff_file.read()
    if gff.startswith("##gff-version"):
        gff_regions = gff.count("sequence-region")
        gff_genes = gff.count("\tgene\t")
        gff_trans = gff.count("\ttranscript\t")
        gff_exons = gff.count("\texon\t")
        gff_cds = gff.count("\tCDS\t")
        log += [f"Supplied GFF file '{input_gff}' was found"]
        log += [
            f"Supplied GFF file contains {gff_regions} sequence regions with:",
            f"    - {gff_genes} genes",
            f"    - {gff_trans} transcripts",
            f"    - {gff_exons} exons",
            f"    - {gff_cds} CDSs",
        ]
    else:
        error += ["Supplied GFF file does not contain a valid '##gff-version' tag"]
    return gff, log, error


if input_database.lower() == "ncbi":
    ncbi_result = getoutput(
        f"datasets summary genome accession {input_assembly} --as-json-lines | "
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
            + f"cd {output_path}; unzip database.zip; rm database.zip"
        )
        copy_command = (
            f"cp {output_path}/ncbi_dataset/data/{refseq_id}/*.fna {output_fasta}; "
            + f"cp {output_path}/ncbi_dataset/data/{refseq_id}/genomic.gff {output_gff}"
        )
        str_out = getoutput(ncbi_command)
        str_cp = getoutput(copy_command)
        # import and check files
        fasta, log, error = check_fasta(output_fasta, log, error)
        gff, log, error = check_gff(output_gff, log, error)

elif input_database.lower() == "manual":
    if not path.exists(input_fasta):
        error += ["The parameter 'fasta' is not a valid path to a FASTA file"]
    elif not path.exists(input_gff):
        error += ["The parameter 'gff' is not a valid path to a GFF/GTF file"]
    else:
        # import and check files
        fasta, log, error = check_fasta(input_fasta, log, error)
        gff, log, error = check_gff(input_gff, log, error)
        # export fasta and gff files
        with open(output_fasta, "w") as fasta_out:
            fasta_out.write(fasta)
        with open(output_gff, "w") as gff_out:
            gff_out.write(gff)
else:
    error += ["The parameter 'database' is none of 'ncbi', 'manual'"]

# print error/log messages
if error:
    print("\n".join(error))
    raise ValueError(
        "Location or format of the supplied genome files was not correct, quitting"
    )
else:
    log += [f"Module finished successfully"]
    log = ["GET_GENOME: " + i for i in log]
    with open(output_log, "w") as log_file:
        log_file.write("\n".join(log))
