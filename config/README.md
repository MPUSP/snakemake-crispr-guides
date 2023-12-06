# Configuration

The workflow can be configured using the `config.yml` file with the following options.

## Parameters

| parameter              | type      | details                                      | default                         |
| ---------------------- | --------- | -------------------------------------------- | ------------------------------- |
| GET_GENOME             |           |                                              |                                 |
| database               | character | one of `ncbi`, `manual`                      | `ncbi`                          |
| assembly               | character | RefSeq ID                                    | `GCF_000006945.2`               |
| fasta                  | path      | optional input                               | `Null`                          |
| gff                    | path      | optional input                               | `Null`                          |
| DESIGN_GUIDES          |           |                                              |                                 |
| target_region          | numeric   | use subset of regions for testing            | `["NC_003277.2"]`               |
| tss_window             | numeric   | upstream/downstream window around TSS        | `[0, 500]`                      |
| circular               | logical   | is the genome circular?                      | `False`                         |
| canonical              | logical   | only canonical PAM sites are included        | `True`                          |
| strands                | character | target `coding`, `template` or `both`        | `both`                          |
| spacer_length          | numeric   | desired length of guides                     | `20`                            |
| guide_aligner          | character | one of `biostrings`, `bowtie`                | `biostrings`                    |
| crispr_enzyme          | character | CRISPR enzyme ID                             | `SpCas9`                        |
| gc_content_range       | numeric   | range of allowed GC content                  | `[30, 70]`                      |
| score_methods          | character | see _crisprScore_ package                    | default scores are listed below |
| score_weights          | numeric   | opt. weights when calculating mean score     | `[1, 1, 1, 1, 1, 1]`            |
| restriction_sites      | character | sequences to omit in entire guide            | `Null`                          |
| bad_seeds              | character | sequences to omit in seed region             | `["ACCCA", "ATACT", "TGGAA"]`   |
| filter_top_n           | numeric   | max number of guides to return               | `10`                            |
| filter_score_threshold | numeric   | mean score to use as lower limit             | `Null`                          |
| filter_multi_targets   | logical   | remove guides that perfectly match >1 target | `True`                          |
| filter_rna             | logical   | remove guides that target e.g. rRNA or tRNA  | `True`                          |
| no_target_controls     | numeric   | number of non-targeting control guides       | `100`                           |
| fiveprime_linker       | character | optionally add 5' linker to each guide       | `Null`                          |
| threeprime_linker      | character | optionally add 3' linker to each guide       | `Null`                          |
| export_as_gff          | logical   | export result table also as `.gff` file      | `False`                         |
| VISUALIZE_GUIDES       |           |                                              |                                 |
| show_examples          | numeric   | number of genes to show guide position       | `10`                            |

## Input data

The workflow requires the following input:

1. An NCBI Refseq ID, e.g. `GCF_000006945.2`. Find your genome assembly and corresponding ID on [NCBI genomes](https://www.ncbi.nlm.nih.gov/data-hub/genome/)
2. OR use a custom pair of `*.fasta` file and `*.gff` file that describe the genome of choice

Important requirements when using custom `*.fasta` and `*.gff` files:

- `*.gff` genome annotation must have the same chromosome/region name as the `*.fasta` file (example: `NC_003197.2`)
- `*.gff` genome annotation must have `gene` and `CDS` type annotation that is automatically parsed to extract transcripts
- `*.gff` genome annotation must have additional qualifiers `Name=...`, `ID=...`, and `Parent=...` for `CDS`s
- all chromosomes/regions in the `*.gff` genome annotation must be present in the `*.fasta` sequence
- but not all sequences in the `*.fasta` file need to have annotated genes in the `*.gff` file

## Further information

Detailed information about exectung the workflow with different command line options can be found in the `README.md`, in the top level dir of this repository.
