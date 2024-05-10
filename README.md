BZ562: Computational Approaches in Molecular Ecology Final Project
==============
## Running the Pariwise Sequentially Markovian Coalescent [(PSMC)](https://github.com/lh3/psmc) on Colorado River Cutthroat Trout (Oncorhynchus clarkii pleuriticus, CRCT) samples. 

- [The Snakeflow](#the-snakeflow)
- [Important files](#important-files)
    - [`sample-info.tsv`](#sample-infotsv)
    - [`config.yaml`](#configyaml)
    - [`Snakefile2`](#snakefile2)
    - [`envs/`](#envs)
    - [`results/`](#results)
- [Major Challenges](#major-challenges)
    - [multiple R1/R2 files per sample](#multiple-r1r2-files-per-sample)
    - [rule mark_duplicates](#rule-mark_duplicates)
    - [Getting PSMC set up](#getting-psmc-set-up)
- [RESULTS!!](#results-1)

## The Snakeflow 
This is a simple Snakemake workflow adapted from Dr. Eric Anderson's [mega-non-model-wgs-snakeflow](https://github.com/eriqande/mega-non-model-wgs-snakeflow).

I processed 4 CRCT samples through this snakeflow for the purpose of the class project. The average read depth for the samples is ~18x (so more like medium coverage WGS data). I aligned these samples to the [Rainbow Trout (O. mykiss) reference genome from NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013265735.2/). 

## Important files

#### `sample-info.tsv`

The 4 samples and their associated read group information can be found in the [sample-info.tsv](https://github.com/noelmclark/crct_snake_proj/blob/main/sample-info.tsv) file. 

#### `config.yaml`

The [config.yaml](https://github.com/noelmclark/crct_snake_proj/blob/main/config.yaml) file contains paths to where my reference fasta and sample-info.tsv are located, as well as a list of all the chormosomes I want to use for future processing, however chromosomes are not needed for the snakeflow to run through the PSMC plot step as I have it right now. 

It also contains the parameters for many of the rules in my Snakefile. 

#### `Snakefile2`

Right now, [Snakefile2](https://github.com/noelmclark/crct_snake_proj/blob/main/Snakefile2) should be used over the older [Snakefile](https://github.com/noelmclark/crct_snake_proj/blob/main/Snakefile). Eventually, I want to separate out the rules into their own files and have the Snakefile include just the paths to the necessary rules as recommended for maximum reproducibility by Snakemake documentation. 

#### `envs/`

This folder contains all of the `.yaml` files for specifying the conda environments including packages and dependencies needed for each rule. 

#### `results/`

This folder includes the PSMC plots that were the ultimate output of the Snakeflow for now... hurrah!!!


## Major Challenges
Here is a copy of my stream of consciousness lab notebook for this project, housed through [Google CoLab](https://colab.research.google.com/drive/1Y0s0Jg9Yp0PwF8AWMJVCd9wf3dBD61Lb#scrollTo=gh-sKP1V9L3a). 

### multiple R1/R2 files per sample

This required adding a unit column to the `sample-info.tsv` file and then using that to identify unique units of different samples until the rule mark_duplicates (GATK MarkDuplicates) stage when you merge them all together. I added the the following code to create a funciton (called get_all_bams_of_common_sample) that creates a list of all the bam files of the same sample (copied from [`common.smk`](https://github.com/eriqande/mega-non-model-wgs-snakeflow/blob/main/workflow/rules/common.smk) file in Eric's [mega-non-model-wgs-snakeflow](https://github.com/eriqande/mega-non-model-wgs-snakeflow).)

```sh
# get all the units of a particular sample
def get_all_bams_of_common_sample(wildcards):
    s=sample_table.loc[(sample_table["sample"] == wildcards.sample)]
    return(expand("results/mapped/{sample}---{unit}.sorted.bam", zip,
        sample = s["sample"].tolist(),
        unit = s["unit"].tolist(),
    ))
```
Additionally, I needed to update the function that is called the expand into the list of the filepaths for each of the R1 and R2 fastq files that so it includes units and deals with some samples that have only one unit per sample. That code is below. 

```sh
def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = sample_table.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}
```
### rule mark_duplicates

At first, the rule mark_duplicates kept failing for the samples with multiple units because the expand rule didn't place commas between file names so it threw a positional error. I played around with the Snakemake wrapper funciton but I had errors with that related to a shell script imbedded in the wrapper so I couldn't really troubleshoot it... 

I went back to the original mark_duplicates rule and it turned out I was missing an ` in the awk command. Once that was working, it was a matter of setting the tmp directory for the program to the right place with enough space and resources for GATK MarkDuplicates to do it's work! 

Here is the final code for the rule mark_duplicates

```sh
rule mark_duplicates:
  input:
    bam=get_all_bams_of_common_sample,
  output:
    bam="results/mkdup/{sample}.bam",
    bai="results/mkdup/{sample}.bai",
    metrics="results/qc/mkdup_metrics/{sample}.metrics"
  conda:
    "envs/gatk.yaml"
  log:
    "results/logs/mark_duplicates/{sample}.log"
  resources:
    mem_mb=112200,
  threads: 30,
  #params:
  #  extra=config["params"]["picard"]["MarkDuplicates"],
  #wrapper:
  #  "v3.9.0/bio/picard/markduplicates"
  shell:
    " BAMS=$(echo {input.bam} | awk '{{for(i=1;i<=NF;i++) printf(\"-I %s \", $i)}}'); "
    " gatk --java-options '-Xmx3740M' MarkDuplicates  "
    "  --CREATE_INDEX --TMP_DIR results/snake-tmp "
    "  $BAMS "
    "  -O {output.bam} "
    "  -M {output.metrics} > {log} 2>&1 "
```

### Getting PSMC set up

PSMC was shockingly easy to set up following the [GitHub documentation](https://github.com/lh3/psmc) for the program. The only error was that the documentation is old and since then the mpileup tool for generating a consensus sequence from a bam file had changed from samtools to bcftoools. Here is the error I got from running it the way the recommend on the PSMC GitHub.

```sh
samtools mpileup -C50 -uf data/genome/OmykA.fasta results/mkdup/C106394.bam | bcftools view -c - | vcfutils.pl vcf2fq -> Activating conda environment: .snakemake/conda/53c9e9ba288a1041ffacd5bee6699548_ mpileup: invalid option -- 'u'

Usage: samtools mpileup [options] in1.bam [in2.bam [...]] Note that using "samtools mpileup" to generate BCF or VCF files has been removed. To output these formats, please use "bcftools mpileup" instead. Error: Could not parse --min-ac - Use of uninitialized value $l in numeric lt (<) at /gpfs/alpine1/scratch/nomclark@colostate.edu/crct_snake_proj/.snakem> Use of uninitialized value $$l in numeric lt (<) at /gpfs/alpine1/scratch/nomclark@colostate.edu/crct_snake_proj/.snakem> Use of uninitialized value $l in numeric lt (<) at /gpfs/alpine1/scratch/nomclark@colostate.edu/crct_snake_proj/.snakem> [Thu May 2 12:47:28 2024] Error in rule psmc_consensus_sequence: jobid: 0 input: results/mkdup/C106394.bam, data/genome/OmykA.fasta output: results/psmc-consensus-sequence/C106394.fq.gz log: results/psmc-consensus-sequence/C106394.log (check log file(s) for error details) conda-env: /gpfs/alpine1/scratch/nomclark@colostate.edu/crct_snake_proj/.snakemake/conda/53c9e9ba288a1041ffacd5bee6> shell: samtools mpileup -C50 -uf data/genome/OmykA.fasta results/mkdup/C106394.bam | bcftools view -c - | vcfutils.pl > (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
```
Below is the code for the rule to generate a psmc consensus sequence needed to run PSMC and plot the results. 

```sh
# rule to get a consensus fastq sequence file for PSMC
# option -d sets and minimum read depth and -D sets the maximum 
# It is recommended to set -d to a third of the average depth and -D to twice
rule psmc_consensus_sequence:
  input:
    bam="results/mkdup/{sample}.bam",
    ref="data/genome/OmykA.fasta",
  output:
    "results/psmc-consensus-sequence/{sample}.fq.gz"
  conda:
    "envs/sambcftools.yaml"
  resources:
    time="23:59:59"
  log:
    "results/psmc-consensus-sequence/{sample}.log"
  shell:
    "bcftools mpileup -C50 -f {input.ref} {input.bam} | bcftools call -c - | " 
    "vcfutils.pl vcf2fq -d 6 -D 36 | gzip > {output} 2> {log}"
```

### RESULTS!!
The rule psmc_plot generated 4 psmc.eps files, one for each sample. I set the generation time `-g` to be 3 years which is standard in the literature for Cutthroat Trout and the mutation rate `-u` to 8.0e-09 based on a [coho salmon PSMC study](https://doi.org/10.1371/journal.pgen.1010918).

```sh
# rule to plot psmc to visualize result
# -u [per-generation mutation rate]
# -g [generation time in years] from https://doi.org/10.1371/journal.pgen.1010918
rule psmc_plot:
  input:
    "results/run-psmc/{sample}.psmc"
  output:
    "results/psmc-plot/{sample}.eps"
  conda:
    "envs/psmc.yaml"
  log:
    "results/psmc-plot/{sample}.log"
  shell:
    "psmc_plot.pl -u 8.0e-09 -g 3 {output} {input}"
```
I took these 4 .eps files into Inkscape to change the colors and make them all fit to the same axes and here is the resulting plot!

![](README_files/all4_psmc.html)<!-- -->