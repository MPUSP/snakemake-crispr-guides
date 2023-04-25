# LOAD PACKAGES
# ------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(devtools)
  library(Biostrings)
  library(GenomicFeatures)
})

devtools::install_github("https://github.com/Bioconductor/BSgenome")
library(BSgenome)

# IMPORT GENOME
# ------------------------------
output <- snakemake@output[["path"]]
messages <- c("importing genome sequence")

genome_fasta <- snakemake@input[["fasta"]]
genome_gff <- snakemake@input[["gff"]]
genome_dna <- Biostrings::readDNAStringSet(genome_fasta)
genome_name <- str_remove_all(names(genome_dna)[1], "^NC\\_[0-9]*\\.[0-9]* |\\,.*")
genome_common <- str_flatten(str_split_1(genome_name, " ")[1:2], " ")
seqinfo_genome <- seqinfo(genome_dna)
isCircular(seqinfo_genome) <- rep_along(seqlevels(seqinfo_genome), FALSE)
genome(seqinfo_genome) <- genome_name
genome_build <- str_remove(read_lines(genome_gff, n_max = 5)[4], "#!genome-build ")

# import genome annotation
txdb <- makeTxDbFromGFF(
  file = genome_gff
)

# check if sequence annotation is identical for sequence and annotation
if (any(seqlevels(seqinfo_genome) != seqlevels(txdb))) {
  seqlevels(seqinfo_genome) <- seqlevels(txdb)
  seqinfo(genome_dna) <- seqinfo_genome
}

# re-import genome annotation with chromosome metadata
txdb <- makeTxDbFromGFF(
  file = genome_gff,
  organism = genome_name,
  chrominfo = seqinfo_genome
)

# create target dirs for output
dir.create(paste0(output, "/seqs_srcdir"))

# create a seed file
txdb_metadata <- metadata(txdb) %>% deframe()
seed_file <- c(
  paste0("Package: BSgenome", str_remove_all(genome_common, " ")),
  paste0("Title: Full genome sequence for ", genome_name),
  paste0("Description: Full genome sequence for ", genome_name),
  "Version: 1.0.0",
  paste0("organism: ", genome_name),
  paste0("common_name: ", genome_common),
  "provider: NCBI",
  paste0("genome: ", "ASM171558v1"),
  paste0("release_date: ", txdb_metadata["Creation time"]),
  paste0("source_url: ", txdb_metadata["Data source"]),
  "organism_biocview: custom_bsgenome",
  "BSgenomeObjname: custom_bsgenome",
  "seqfile_name: single_sequences.2bit",
  # paste0("seqnames: ", "'single_sequences'"),#c('", paste(txdb$user_seqlevels, collapse = "', '"), "')"),
  paste0("seqs_srcdir: ", output, "/seqs_srcdir")
)

# export seed file
write_lines(
  file = paste0(output, "/bsgenome_seed_file.txt"),
  x = seed_file
)

# export sequence file
for (seqs in seqlevels(genome_dna)) {
  genome_dna[seqs] %>%
    writeXStringSet(paste0(output, "/seqs_srcdir/", seqs, ".fa"), format = "fasta")
}

forgeSeqFiles(
  provider = "NCBI",
  genome = "ASM171558v1",
  seqnames = txdb$user_seqlevels,
  suffix = ".fa",
  seqs_srcdir = paste0(output, "/seqs_srcdir"),
  seqs_destdir = paste0(output, "/seqs_srcdir"),
  ondisk_seq_format = "2bit",
  verbose = TRUE
)

forgeBSgenomeDataPkg(
  x = paste0(output, "/bsgenome_seed_file.txt"),
  destdir = output
)

# # install BSgenome package
# devtools::install(paste0(
#   output, "/package/BSgenome",
#   str_remove_all(genome_common, " ")
# ))


# # export log file
# write_lines(
#   file = snakemake@log[["path"]],
#   x = paste0("DESIGN_GUIDES: ", messages)
# )
