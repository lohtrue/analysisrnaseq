# lohtrue/analysisrnaseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [2026-10-23]

Initial release of lohtrue/analysisrnaseq, created with the [nf-core](https://nf-co.re/) template (nf-core/tools 4.1.0).

### `Added`

- Samplesheet as input with single- and paired-end FASTQ files and strandedness (`forward`, `reverse`, `unstranded`), validated with nf-schema
- Optional subsampling of reads with `--subsample` (fq subsample, fixed seed)
- Read quality control with FastQC
- Adapter/quality trimming with fastp
- Decompression of reference files with gunzip - if gzipped
- Alignment to the reference genome with HISAT2, including extraction of splice sites and index building 
- Sorting of alignments with samtools sort
- Duplicate marking with samtools markdup
- Alignment statistics with samtools flagstat
- Gene-level read counting with featureCounts
- Local module `TPM_TABLE` creates TPM expression table (genes × samples)
- Local module `PLOT_TPM` plots of the top 20 genes, the expression distribution and the sample correlation, included in the MultiQC report
- MultiQC report summarising all QC results and software versions
- Resource limit parameters `--max_cpus`, `--max_memory` and `--max_time` (defaults: 4 CPUs, 6 GB, 6 h) for running on laptops
- Test profiles `test` (mouse chr19) and `test_full` (yeast data from nf-core/test-datasets)
- Documentation: README with metro map, usage and output documentation, citations

### `Dependencies`

| Tool          | Version |
| ------------- | ------- |
| fq            | 0.12.0  |
| FastQC        | 0.12.1  |
| fastp         | 1.3.6   |
| HISAT2        | 2.2.3   |
| SAMtools      | 1.24    |
| Subread       | 2.1.1   |
| MultiQC       | 1.35    |
| matplotlib    | 3.11.2  |