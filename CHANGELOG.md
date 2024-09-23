# Changelog

## [1.3.0](https://github.com/MPUSP/snakemake-crispr-guides/compare/v1.2.0...v1.3.0) (2024-09-23)


### Features

* changed license to MIT, closes [#30](https://github.com/MPUSP/snakemake-crispr-guides/issues/30) ([c1627fb](https://github.com/MPUSP/snakemake-crispr-guides/commit/c1627fb059252a07b2ce506a50be910e202c3dae))
* minor updates to README ([516a728](https://github.com/MPUSP/snakemake-crispr-guides/commit/516a728c70af8e7c001851d6e10fafb85b00be88))


### Bug Fixes

* bug when restriction sites is NULL ([d26fc02](https://github.com/MPUSP/snakemake-crispr-guides/commit/d26fc02b8384528373aad866434571225ea02512))
* bump container version for next release ([d893c6b](https://github.com/MPUSP/snakemake-crispr-guides/commit/d893c6be7f7315646b835139756d68d5c58d9bba))
* handling restriction sites, closes [#36](https://github.com/MPUSP/snakemake-crispr-guides/issues/36) ([675747d](https://github.com/MPUSP/snakemake-crispr-guides/commit/675747d76d980f01514fb8d420dfb3312a12ca29))
* minor issues, closes [#34](https://github.com/MPUSP/snakemake-crispr-guides/issues/34) ([8b7b5bb](https://github.com/MPUSP/snakemake-crispr-guides/commit/8b7b5bbf2da92c1733cf77eb35a90b68c5344a52))
* removed conda default channel from envs, closes [#38](https://github.com/MPUSP/snakemake-crispr-guides/issues/38) ([0362a0b](https://github.com/MPUSP/snakemake-crispr-guides/commit/0362a0bea66e1e38bdd698f489e2f7abdf27df70))

## [1.2.0](https://github.com/MPUSP/snakemake-crispr-guides/compare/v1.1.0...v1.2.0) (2024-02-21)


### Features

* added possibility to run with singularity container ([9d4598e](https://github.com/MPUSP/snakemake-crispr-guides/commit/9d4598ea1fb21ff262250d961da8efd23de96a46))
* option to add padding between CDS and intergenic region, closes [#25](https://github.com/MPUSP/snakemake-crispr-guides/issues/25) ([de3f8bf](https://github.com/MPUSP/snakemake-crispr-guides/commit/de3f8bfd551b4bd960648baed3e86cc87235b7fd))


### Bug Fixes

* added .vscode ([d4a15e0](https://github.com/MPUSP/snakemake-crispr-guides/commit/d4a15e0698e93c608831fefd7de6cd650736591d))
* bug in export of NTC guides when all are filtered out ([bf97758](https://github.com/MPUSP/snakemake-crispr-guides/commit/bf97758a440ac7c31db21ecd8678ad2289e160c8))
* overhauled plotting of intergenic regions, closes [#28](https://github.com/MPUSP/snakemake-crispr-guides/issues/28) ([da7135c](https://github.com/MPUSP/snakemake-crispr-guides/commit/da7135cc0d92376e69047343280106a388fb9d44))
* updated container, bumped version ([cc9dcf9](https://github.com/MPUSP/snakemake-crispr-guides/commit/cc9dcf92f3201cb3d52cb4a08d969e5a1ab591d8))

## [1.1.0](https://github.com/MPUSP/snakemake-crispr-guides/compare/v1.0.0...v1.1.0) (2024-01-12)


### Features

* added config README and updated workflow catalog options ([13649a1](https://github.com/MPUSP/snakemake-crispr-guides/commit/13649a179c621411a5df830d2d7a9c35626a25ca))
* added image with pipeline overview ([69a816c](https://github.com/MPUSP/snakemake-crispr-guides/commit/69a816c7bc07c8db3a1f479c49526d71a10bcc64))
* added option to export guide list as gff file, closes [#16](https://github.com/MPUSP/snakemake-crispr-guides/issues/16) ([470025c](https://github.com/MPUSP/snakemake-crispr-guides/commit/470025c666c193cbe165bb39f6de5964391c2025))
* major and changes to support tiling based guide design ([478c68f](https://github.com/MPUSP/snakemake-crispr-guides/commit/478c68f87673baeb533fcbd2dfb184267967a492))
* more flexibility with chromosome names, updated envs, closes [#15](https://github.com/MPUSP/snakemake-crispr-guides/issues/15) ([edd1d16](https://github.com/MPUSP/snakemake-crispr-guides/commit/edd1d1601688372d7f84bf3bb3866ed96377e442))
* parsing of GFF file improved considerably by using bcbio GFF lib, [#15](https://github.com/MPUSP/snakemake-crispr-guides/issues/15) ([80c66ff](https://github.com/MPUSP/snakemake-crispr-guides/commit/80c66ff35eabda23d190609764332477e4e54fef))
* smarter retrieveal of genome name and tax ID; closes [#19](https://github.com/MPUSP/snakemake-crispr-guides/issues/19) ([afdc2e3](https://github.com/MPUSP/snakemake-crispr-guides/commit/afdc2e3827ee46fce77eca0c838f4a5cc9d746af))


### Bug Fixes

* added --cores to lint command ([c9fd2a5](https://github.com/MPUSP/snakemake-crispr-guides/commit/c9fd2a589f7eee0ae23943839f38bb0272899005))
* added missing example pics ([7a935a8](https://github.com/MPUSP/snakemake-crispr-guides/commit/7a935a8a4012618ec827444f7436603bb5f05002))
* added option to specify GFF sources in config file; closes [#21](https://github.com/MPUSP/snakemake-crispr-guides/issues/21) ([b02e869](https://github.com/MPUSP/snakemake-crispr-guides/commit/b02e869945abdee5b2c96f2c045862fb7468768a))
* availability of tx_list when target is empty ([2e0a248](https://github.com/MPUSP/snakemake-crispr-guides/commit/2e0a248e852810440445088e70c21421cf5e92d8))
* bumped versions and removed unnecessary dependency ([2a1975a](https://github.com/MPUSP/snakemake-crispr-guides/commit/2a1975abf0e233d39aa2c7a40c9c2dad46edf42e))
* extended documentation for 'usage' page on sm-workflow-catalog ([bd6b4fe](https://github.com/MPUSP/snakemake-crispr-guides/commit/bd6b4fef6cd3ea4c1520ed4cbde4b8a337e14d0c))
* removed 1 branch from CI test workflow ([f3fe343](https://github.com/MPUSP/snakemake-crispr-guides/commit/f3fe343bee4ccaafa1841e0b8af4e92e5438916e))
* removed a comment ([ac5a6e2](https://github.com/MPUSP/snakemake-crispr-guides/commit/ac5a6e2bfb8ffabede32caaa9f2dd98c14ce035f))
* removed additional flag fields entirely ([53319a0](https://github.com/MPUSP/snakemake-crispr-guides/commit/53319a0c279ed29f322372d216a9f230b57a33b4))
* restored original empty fields for flags ([01ef026](https://github.com/MPUSP/snakemake-crispr-guides/commit/01ef0265e24d22361895fc7edbe21497d67e1d69))

## 1.0.0 (2023-07-13)


### Features

* add major module to find, score and filter guide RNAs ([4f93545](https://github.com/MPUSP/snakemake-crispr-guides/commit/4f935454f1496fef2872fe26771fcf959c4fe5ff))
* added basic R markdown report module ([3950bcd](https://github.com/MPUSP/snakemake-crispr-guides/commit/3950bcde6a23ff600ef0c60d2301567087e2c6cf))
* added CI test dir ([bdfe96c](https://github.com/MPUSP/snakemake-crispr-guides/commit/bdfe96c8400cf13edd80e91b1888328d8cc7db38))
* added dev branch to GH actions ([303de61](https://github.com/MPUSP/snakemake-crispr-guides/commit/303de61ca37ad63d3a39de119096df10206ea94f))
* added function to exclude RNAs from guide design ([63f0b0e](https://github.com/MPUSP/snakemake-crispr-guides/commit/63f0b0ec2ce995c8d23ebf2489516cae63685425))
* added functionality to design random control guides ([a32185b](https://github.com/MPUSP/snakemake-crispr-guides/commit/a32185ba426b6de48c4c5c249620ef1e018a1f8f))
* added functionality to handle and visualize guides with multiple targets. Closes [#11](https://github.com/MPUSP/snakemake-crispr-guides/issues/11) ([ac4cd79](https://github.com/MPUSP/snakemake-crispr-guides/commit/ac4cd79088a461be873be80d59ddba017115ce5c))
* added log entry with number of non-targeted genes ([1b1653f](https://github.com/MPUSP/snakemake-crispr-guides/commit/1b1653f2c464258823851487ec32b7b0b13e476c))
* added module to convert html report to pdf ([37f30a4](https://github.com/MPUSP/snakemake-crispr-guides/commit/37f30a4a35a55875d9bef04b308b7f0c111fc359))
* added module to generate BSgenome object on the fly, for use with bowtie, closes [#4](https://github.com/MPUSP/snakemake-crispr-guides/issues/4) ([c344d0d](https://github.com/MPUSP/snakemake-crispr-guides/commit/c344d0d70ea42b147b76df41d8dbe99c91526ab1))
* added module to make bowtie2 index ([cf11e73](https://github.com/MPUSP/snakemake-crispr-guides/commit/cf11e73821fc56a4d30f12e65f43d521686458a1))
* added more TSS features to result table, and new viz figures ([8ffa6a2](https://github.com/MPUSP/snakemake-crispr-guides/commit/8ffa6a24fa74ba593a40f36222850a7e0f994b48))
* added option to add arbitrary linkers to guide RNA ends ([d931cb7](https://github.com/MPUSP/snakemake-crispr-guides/commit/d931cb722c00cb7e75f0a673f46a400cf69dbbb7))
* added optional weights to calculate mean score ([5678c38](https://github.com/MPUSP/snakemake-crispr-guides/commit/5678c38689dcd374aa483074f2f191d3dcea22e2))
* added output of failed guide RNA design to report ([4019861](https://github.com/MPUSP/snakemake-crispr-guides/commit/4019861fecd5a3b5b9ea9dccb953261ce6a7e109))
* added percentages to report plots ([da8b5d9](https://github.com/MPUSP/snakemake-crispr-guides/commit/da8b5d9fdb8ef8e8b5185297bece70b6e246881b))
* added score for G eenrichment, closes [#5](https://github.com/MPUSP/snakemake-crispr-guides/issues/5) ([a1b35a1](https://github.com/MPUSP/snakemake-crispr-guides/commit/a1b35a136ffc79fcf2c979e8ec1f784d00602e5d))
* changed license to GPL ([9ce508a](https://github.com/MPUSP/snakemake-crispr-guides/commit/9ce508a30f8b9ba3ee00fcd5d36a298fe9ab58d2))
* export table of targets where guide design failed ([e377f4b](https://github.com/MPUSP/snakemake-crispr-guides/commit/e377f4b8a29ca88803c216c3d8a57fbf0d34b997))
* filter for guides that fall within TSS window _and_ transcript region + visualization ([39e9cea](https://github.com/MPUSP/snakemake-crispr-guides/commit/39e9ceadebb3211c879a2c7ccfeac664ea3befb3))
* finished module for obtaining genomes from NCBI or user ([5c1fa2b](https://github.com/MPUSP/snakemake-crispr-guides/commit/5c1fa2b97c92e9aa362d05a33d117f686ae50e54))
* first draft of module to generate bsgenome object on the fly ([a5f8b77](https://github.com/MPUSP/snakemake-crispr-guides/commit/a5f8b77633d6d0ff0dc4751563a2eb65e614eb6e))
* fixed small bugs to improve usage of custom fasta/gff files ([ab76a25](https://github.com/MPUSP/snakemake-crispr-guides/commit/ab76a25b153d70415faae4113e8bf1044a3662af))
* improved prediction around TSS and added TSS dist score, closes [#3](https://github.com/MPUSP/snakemake-crispr-guides/issues/3) ([290b675](https://github.com/MPUSP/snakemake-crispr-guides/commit/290b67569d1702da1068bca6cdeea5396b95f98c))
* included visualization of guide RNA location, and option for strand preference ([d345099](https://github.com/MPUSP/snakemake-crispr-guides/commit/d3450996b001d579078bf10885adad5a41aa9705))
* simplified handling of result files ([ddcb925](https://github.com/MPUSP/snakemake-crispr-guides/commit/ddcb925fe87a4ab2d1dd6f9c911543a14d44b44f))


### Bug Fixes

* added missing documentation for paramter filter_multi_targets ([32476e5](https://github.com/MPUSP/snakemake-crispr-guides/commit/32476e57be7fba97112bd29aa946adfff8520395))
* added r-base environment to conda list ([e2d7cbf](https://github.com/MPUSP/snakemake-crispr-guides/commit/e2d7cbfa92e6e6d884696afec41addf06fc88140))
* added variable TSS windows, closes [#9](https://github.com/MPUSP/snakemake-crispr-guides/issues/9) ([014e6ec](https://github.com/MPUSP/snakemake-crispr-guides/commit/014e6ec7ca24071b415ad1f8412d081ebfe0a447))
* changed syntax for env defs, use dry run for test ([32f5e82](https://github.com/MPUSP/snakemake-crispr-guides/commit/32f5e82525e514e455f66aef948b26e63d405db2))
* conda env error nothing provides libgcc-ng ([ffc9cec](https://github.com/MPUSP/snakemake-crispr-guides/commit/ffc9cec6129459b40b8672a05bd48f1f26a9f0e6))
* fixed small syntax error for env ([bd1775b](https://github.com/MPUSP/snakemake-crispr-guides/commit/bd1775baec6880a7ee31b371a0612d7b853b1ba9))
* improvements to score correlation plot ([7afc9cb](https://github.com/MPUSP/snakemake-crispr-guides/commit/7afc9cb488f3b053766b7739ee4d69dfeef71927))
* removed typo in dry snakemake test command ([8e3e5c0](https://github.com/MPUSP/snakemake-crispr-guides/commit/8e3e5c0b94769a139710ad11171a7448d440d67b))
* replaced bowtie2 with bowtie1 as a req for crisprDesign module ([f8cd5fb](https://github.com/MPUSP/snakemake-crispr-guides/commit/f8cd5fbd985ab2a4ed80c73c66fdbf3157c20567))
* small improvements to report and docs ([d24c36c](https://github.com/MPUSP/snakemake-crispr-guides/commit/d24c36c7889efc4dba083a75c8ea859525d40a29))
* small update to guide RNA drawing ([18eac19](https://github.com/MPUSP/snakemake-crispr-guides/commit/18eac19aa440b3c8b6bf89ac9bf018d2a892c16f))
* some corrections to pipeline overview ([2e519ed](https://github.com/MPUSP/snakemake-crispr-guides/commit/2e519ed8a1b7b98b4b45cbdef4890128dbeb0ae2))
