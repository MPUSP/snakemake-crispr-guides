# Changelog

## 1.0.0 (2023-04-26)


### Features

* add major module to find, score and filter guide RNAs ([4f93545](https://github.com/MPUSP/snakemake-crispr-guides/commit/4f935454f1496fef2872fe26771fcf959c4fe5ff))
* added basic R markdown report module ([3950bcd](https://github.com/MPUSP/snakemake-crispr-guides/commit/3950bcde6a23ff600ef0c60d2301567087e2c6cf))
* added CI test dir ([bdfe96c](https://github.com/MPUSP/snakemake-crispr-guides/commit/bdfe96c8400cf13edd80e91b1888328d8cc7db38))
* added dev branch to GH actions ([303de61](https://github.com/MPUSP/snakemake-crispr-guides/commit/303de61ca37ad63d3a39de119096df10206ea94f))
* added module to convert html report to pdf ([37f30a4](https://github.com/MPUSP/snakemake-crispr-guides/commit/37f30a4a35a55875d9bef04b308b7f0c111fc359))
* added module to make bowtie2 index ([cf11e73](https://github.com/MPUSP/snakemake-crispr-guides/commit/cf11e73821fc56a4d30f12e65f43d521686458a1))
* added more TSS features to result table, and new viz figures ([8ffa6a2](https://github.com/MPUSP/snakemake-crispr-guides/commit/8ffa6a24fa74ba593a40f36222850a7e0f994b48))
* added optional weights to calculate mean score ([5678c38](https://github.com/MPUSP/snakemake-crispr-guides/commit/5678c38689dcd374aa483074f2f191d3dcea22e2))
* added score for G eenrichment, closes [#5](https://github.com/MPUSP/snakemake-crispr-guides/issues/5) ([a1b35a1](https://github.com/MPUSP/snakemake-crispr-guides/commit/a1b35a136ffc79fcf2c979e8ec1f784d00602e5d))
* changed license to GPL ([9ce508a](https://github.com/MPUSP/snakemake-crispr-guides/commit/9ce508a30f8b9ba3ee00fcd5d36a298fe9ab58d2))
* filter for guides that fall within TSS window _and_ transcript region + visualization ([39e9cea](https://github.com/MPUSP/snakemake-crispr-guides/commit/39e9ceadebb3211c879a2c7ccfeac664ea3befb3))
* finished module for obtaining genomes from NCBI or user ([5c1fa2b](https://github.com/MPUSP/snakemake-crispr-guides/commit/5c1fa2b97c92e9aa362d05a33d117f686ae50e54))
* fixed small bugs to improve usage of custom fasta/gff files ([ab76a25](https://github.com/MPUSP/snakemake-crispr-guides/commit/ab76a25b153d70415faae4113e8bf1044a3662af))
* improved prediction around TSS and added TSS dist score, closes [#3](https://github.com/MPUSP/snakemake-crispr-guides/issues/3) ([290b675](https://github.com/MPUSP/snakemake-crispr-guides/commit/290b67569d1702da1068bca6cdeea5396b95f98c))
* included visualization of guide RNA location, and option for strand preference ([d345099](https://github.com/MPUSP/snakemake-crispr-guides/commit/d3450996b001d579078bf10885adad5a41aa9705))
* simplified handling of result files ([ddcb925](https://github.com/MPUSP/snakemake-crispr-guides/commit/ddcb925fe87a4ab2d1dd6f9c911543a14d44b44f))
