## Running the workflow

### Input data

The workflow requires the following input:

1. An NCBI Refseq ID, e.g. `GCF_000006945.2`. Find your genome assembly and corresponding ID on [NCBI genomes](https://www.ncbi.nlm.nih.gov/data-hub/genome/)
2. OR use a custom pair of `*.fasta` file and `*.gff` file that describe the genome of choice

Important requirements when using custom `*.fasta` and `*.gff` files:

- `*.gff` genome annotation must have the same chromosome/region name as the `*.fasta` file (example: `NC_003197.2`)
- `*.gff` genome annotation must have `gene` and `CDS` type annotation that is automatically parsed to extract transcripts
- `*.gff` genome annotation must have additional qualifiers `Name=...`, `ID=...`, and `Parent=...` for `CDS`s
- all chromosomes/regions in the `*.gff` genome annotation must be present in the `*.fasta` sequence
- but not all sequences in the `*.fasta` file need to have annotated genes in the `*.gff` file

### Execution

To run the workflow from command line, change the working directory.

```
cd /path/to/snakemake-crispr-guides
```

Adjust the global and module-specific options in the default config file `config/config.yml`.
Before running the entire workflow, you can perform a dry run using:

```
snakemake --dry-run
```

To run the complete workflow with test files, execute the following command. The definition of the number of compute cores is mandatory.

```
snakemake --cores 10 --use-conda
```

To supply a custom config file and/or use options that override the defaults, run the workflow like this:

```
snakemake --cores 10 --use-conda \
  --configfile 'config/my_config.yml' \
  --config \
  option='input'
```

### Parameters

This table lists all parameters that can be used to run the workflow.

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

### Off-target scores

The pipeline maps each guide RNA to the target genome and -by default- counts the number of alternative alignments with 1, 2, 3, or 4 mismatches. All guide RNAs that map to any other position including up to 4 allowed mismatches are removed.
An exception to this rule is made for guides that perfectly match multiple targets when the `filter_multi_targets` is set to `False` (default: `True`). The reasoning behind this rule is that genomes often contain duplicated genes/targets, and the default but sometimes undesired behavior is to remove all guides targeting the two or more duplicates. If set to `False`, these guides will not be removed and duplicated genes will be targeted even if they are located at different sites. 

### On-target scores

The list of available on-target scores in the [R crisprScore package](https://github.com/crisprVerse/crisprScore) is larger than the different scores included by default. It is important to note that the computation of some scores does not necessarily make sense for the design of every CRISPR library. For example, several scores were obtained from analysis of Cas9 cutting efficiency in human cell lines. For such scores it is questionable if they are useful for the design of a different type of library, for example a dCas9 CRISPR inhibition library for bacteria.

Another good reason to exclude some scores are the computational resources they require. Particularly deep learning-derived scores are calculated by machine learning models that require both a lot of extra resources in terms of disk space (downloaded and installed *via* `basilisk` and `conda` environments) and processing power (orders of magnitude longer computation time).

Users can look up all available scores on the [R crisprScore github page](https://github.com/crisprVerse/crisprScore) and decide which ones should be included. In addition, the default behavior of the pipeline is to compute an average score and select the top N guides based on it. The average score is the _weighted mean_ of all single scores and the `score_weights` can be defined in the `config/config.yml` file. If a score should be excluded from the ranking, it's weight can simply be set to zero.

The default scores are:

- `ruleset1`, `ruleset3`, `crisprater`, and `crisprscan` from the `crisprScore` package
- `tssdist` as an additional score representing the relative distance to the promoter. Only relevant for CRISRPi repression
- `genrich` as an additional score representing the `G` enrichment in the -4 to -14 nt region of a spacer ([Miao & Jahn et al., 2023](https://www.biorxiv.org/content/10.1101/2023.02.13.528328v1)). Only relevant for CRISPRi repression

### Strand specificity

The strand specificity is important for some CRISPR applications. In contrast to the `crisprDesign` package, functions were added to allow the design of guide RNAs that target either both strands, or just the coding (non-template) strand, or the template strand. This can be defined with the `strands` parameter in the config file.

- For CRISPRi (inhibition) experiments, the literature recommends to target the **coding strand for the CDS** or **both strands for the promoter** ([Larson et al., Nat Prot, 2013](http://dx.doi.org/10.1038/nprot.2013.132))
- this pipeline will automatically filter guides for the chosen strand
- for example, if only guides for the coding (non-template) strand are desired, genes on the "+" strand will be targeted with reverse-complement guides ("-"), and genes on the "-" strand with "+" guides.

### Random control guides

The pipeline includes the option to design random control guides. These are simply random nucleotide sequences with the same length as specified for the actual guide RNAs. The control guides are named `NTC_<number>` and exported in a separate table, by default as `results/design_guides/guideRNAs_ntc.csv`. Some very reduced checks are done for these guides, such as off-target binding. However most on-target checks are omitted for these guides as they have no defined binding site, strand, or other typical guide properties. Linkers are added when specified.
