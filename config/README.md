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

### Target type

One of the most important options is to specify the _type of target_ with the `target_type` parameter. The pipeline can generate up to three different types of guide RNAs:

1. guides for **targets** - these are typically genes, promoters or other annotated genetic elements determined from the supplied GFF file. The pipeline will try to find the best guides by position and score targeting the defined window around the start of the gene/feature (parameter `tss_window`). The number of guides is specified with `filter_best_per_gene`.
2. guides for **intergenic regions** - for non-annotated regions (or in the absence of any targets), the pipeline attempts to design guide RNAs using a 'tiling' approach. This means that the supplied genome is subdivided into 'tiles' (bins) of width `tiling_window`, and the best guide RNAs per window are selected. The number of guides is specified with `filter_best_per_tile`.
3. guides **not targeting anything** - this type of guide RNAs is most useful as negative control, in order to gauge the effect of the genetic background on mutant selection without targeting a gene. These guides are random nucleotide sequences with the same length as the target guide RNAs. The no-target control guides are named `NTC_<number>` and exported in a separate table (`results/filter_guides/guideRNAs_ntc.csv`). Some very reduced checks are done for these guides, such as off-target binding. mMst on-target checks are omitted for these guides as they have no defined binding site, strand, or other typical guide properties.

### Off-target scores

The pipeline maps each guide RNA to the target genome and -by default- counts the number of alternative alignments with 1, 2, 3, or 4 mismatches. All guide RNAs that map to any other position including up to 4 allowed mismatches are removed.
An exception to this rule is made for guides that perfectly match multiple targets when the `filter_multi_targets` is set to `False` (default: `True`). The reasoning behind this rule is that genomes often contain duplicated genes/targets, and the default but sometimes undesired behavior is to remove all guides targeting the two or more duplicates. If set to `False`, these guides will not be removed and duplicated genes will be targeted even if they are located at different sites.

### On-target scores

The list of available on-target scores in the [R crisprScore package](https://github.com/crisprVerse/crisprScore) is larger than the different scores included by default. It is important to note that the computation of some scores does not necessarily make sense for the design of every CRISPR library. For example, several scores were obtained from analysis of Cas9 cutting efficiency in human cell lines. For such scores it is questionable if they are useful for the design of a different type of library, for example a dCas9 CRISPR inhibition library for bacteria.

Another good reason to exclude some scores are the computational resources they require. Particularly deep learning-derived scores are calculated by machine learning models that require both a lot of extra resources in terms of disk space (downloaded and installed _via_ `basilisk` and `conda` environments) and processing power (orders of magnitude longer computation time).

Users can look up all available scores on the [R crisprScore GitHub page](https://github.com/crisprVerse/crisprScore) and decide which ones should be included. In addition, the default behavior of the pipeline is to compute an average score and select the top N guides based on it. The average score is the _weighted mean_ of all single scores and the `score_weights` can be defined in the `config/config.yml` file. If a score should be excluded from the ranking, it's weight can simply be set to zero.

The default scores are:

- `ruleset1`, `ruleset3`, `crisprater`, and `crisprscan` from the `crisprScore` package
- `tssdist` as an additional score representing the relative distance to the promoter. Only relevant for CRISRPi repression
- `genrich` as an additional score representing the `G` enrichment in the -4 to -14 nt region of a spacer ([Miao & Jahn et al., 2023](https://www.biorxiv.org/content/10.1101/2023.02.13.528328v1)). Only relevant for CRISPRi repression

### Strand specificity

The strand specificity is important for some CRISPR applications. In contrast to the `crisprDesign` package, functions were added to allow the design of guide RNAs that target either both strands, or just the coding (non-template) strand, or the template strand. This can be defined with the `strands` parameter in the config file.

- For CRISPRi (inhibition) experiments, the literature recommends to target the **coding strand for the CDS** or **both strands for the promoter** ([Larson et al., Nat Prot, 2013](http://dx.doi.org/10.1038/nprot.2013.132))
- this pipeline will automatically filter guides for the chosen strand
- for example, if only guides for the coding (non-template) strand are desired, genes on the "+" strand will be targeted with reverse-complement guides ("-"), and genes on the "-" strand with "+" guides.
