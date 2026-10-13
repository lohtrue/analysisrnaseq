# lohtrue/analysisrnaseq: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

lohtrue/analysisrnaseq is a lightweight RNA-analysis pipeline that turns sequencing reads into a table of gene expression values. For each sample of the samplesheet, every read gets quality-checked, trimmed and aligned to a reference genome. The reads per gene get counted after marking duplicate reads. Finally, the gene counts are converted to TPM(transcripts per million) and summarized in a table, as well as visualized with overview plots and a MultiQC report summarizes the QC-results.

The pipeline is designed for educational purposes with the goal of easy execution on laptops with limited resources. To run lohtrue/analysisrnaseq Nextflow, Docker, a samplesheet listing the FASTQ-files, a reference genome (FASTA) and an annotation file (GTF) are needed.

## Samplesheet input

Before running the pipeline, a samplesheets with information about the samples (FASTQ-files) needs to be created. It has to be a comma-separated file with 4 columns and a header row as shown in the examples below. To specify the location of the samplesheet the parameter --input must be used:

```bash
--input '[path to samplesheet file]'
```
The samplesheet must contain the following 4 columns: 

| Column         | Description |
| -------------- | ----------- |
| `sample`       | Sample name without spaces. Used as the sample ID in all output files and as the column name in the TPM table. |
| `fastq_1`      | Path to the (gzipped) FASTQ file with read 1 (`.fastq.gz` or `.fq.gz`). |
| `fastq_2`      | Path to the (gzipped) FASTQ file with read 2. Leave empty for single-end data. |
| `strandedness` | Strandedness: `forward`, `reverse` or `unstranded`. Needed for correct alignmentand counting . |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

To check that everything is set up correctly, you can first run the pipeline on the small test dataset: 

```bash
nextflow run lohtrue/analysisrnaseq -r dev -profile test,docker --outdir results_test
```

A typical command for your own data looks like this:

```bash
nextflow run lohtrue/analysisrnaseq -r dev \
    -profile docker \
    --input samplesheet.csv \
    --fasta genome.fa.gz \
    --gtf annotation.gtf.gz \
    --outdir results
```

The pipeline creates the following files in your working directory:

```bash
work            # Intermediate files of all tasks (can be deleted after the run) 
<OUTDIR>        # Final results (defined with --outdir)
.nextflow.log   # Log file from Nextflow
```

### Using a params file

Instead of writing all parameters into the command, they can be stored in a YAML file and passed with `-params-file`. 

```bash
nextflow run lohtrue/analysisrnaseq -r dev -profile docker -params-file params.yaml
```

```yaml title="params.yaml"
input: "samplesheet.csv"
fasta: "genome.fa.gz"
gtf: "annotation.gtf.gz"
outdir: "results"
subsample: 0.1
max_memory: "6.GB"
```

> [!WARNING]
> Do not use `-c <file>` to set parameters. Config files given with `-c` should only be used for process settings such as tool arguments.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull lohtrue/analysisrnaseq -r dev
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline.

Use `-r` (one hyphen) to choose the version: either a branch (e.g. `-r dev`) or, once releases are available, a version tag (e.g. `-r 1.0.0`). The version used is logged in the reports, for example at the bottom of the MultiQC report, so you can see later which version produced your results.

To further assist in reproducibility, you can share and reuse [parameter files](#using-a-params-file) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).

Only `docker` has been tested with this pipeline. Conda is not supported, because the local modules (`tpm_table`, `plot_tpm`) do not provide a Conda environment.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource limits

By default, the pipeline is set up for computers with limited resources; each task uses at most 4 CPUs, 6 GB of memory and 6 hours. These limits can be changed with the following parameters:

| Parameter      | Default | Description                     |
| -------------- | ------- | ------------------------------- |
| `--max_cpus`   | `4`     | Maximum number of CPUs per task |
| `--max_memory` | `6.GB`  | Maximum memory per task         |
| `--max_time`   | `6.h`   | Maximum runtime per task        |


### Subsampling

With `--subsample`, only a specified fraction of the reads of each sample is used, e.g. `--subsample 0.1` keeps 10% of the reads. The value must be between 0 and 1. A fixed seed is used, so the same reads are selected in every run. By default, all reads are used.

### Example: running on laptop

The following command runs the pipeline on a laptop with 2 CPUs and 5 GB RAM, using only 10 % of the reads of each sample: 

```bash
nextflow run lohtrue/analysisrnaseq -r dev \
    -profile docker \
    --input samplesheet.csv \
    --fasta genome.fa.gz \
    --gtf annotation.gtf.gz \
    --outdir results \
    --subsample 0.1 \
    --max_cpus 2 \
    --max_memory 5.GB \
    --max_time 12.h
```

### Custom tool arguments

The arguments passed to each tool are set with `ext.args` in [`conf/modules.config`](../conf/modules.config). To change them without editing the pipeline, create your own config file and pass it with `-c`:

```groovy title="custom.config"
process {
    withName: 'FASTP' {
        ext.args = '--qualified_quality_phred 20'
    }
}
```

## Running in the background

Nextflow must keep running until the pipeline has finished. With `-bg`, Nextflow runs in the background and writes its output to `.nextflow.log`. On a laptop, make sure it does not go to sleep during the run.

## Nextflow memory requirements

Nextflow itself runs in Java and can use a lot of memory. To limit this, add the following line to your shell configuration (`~/.zshrc` on macOS, `~/.bashrc` on Linux):

```bash
export NXF_OPTS='-Xms1g -Xmx2g'
```