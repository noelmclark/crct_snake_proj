# start with some stuff to get the Snakemake version
from snakemake import __version__ as snakemake_version
smk_major_version = int(snakemake_version[0])



# import modules as needed. Also import the snakemake
# internal command to load a config file into a dict
import pandas as pd

if smk_major_version >= 8:
  from snakemake.common.configfile import _load_configfile
else:
  from snakemake.io import _load_configfile



# define rules that don't need to be run on a compute node.
# i.e. those that can be run locally.
localrules: all, genome_faidx, genome_dict



### Get a dict named config from crct_snake_proj/config.yaml
configfile: "config.yaml"




### Get the sample info table read into a pandas data frame
sample_table=pd.read_table(config["sample_info"], dtype="str").set_index(
    "sample", drop=False
)



### Transfer values from the yaml and tabular config to
### our familiar lists, SAMPLES and CHROMOS
# Populate our SAMPLES list from the sample_table using a little
# pandas syntax
SAMPLES=sample_table["sample"].unique().tolist()

# Create a new SAMPLESID list from the sample_table using the same
# pandas syntax as for SAMPLES
SAMPLEID=sample_table["sampleid"].unique().tolist()

# Define CHROMOS from the values in the config file
CHROMOS=config["chromos"]



### Input Functins that use the tabular sample_info
# define a function to get the fastq path from the sample_table. This
# returns it as a dict, so we need to unpack it in the rule
def get_fastqs(wildcards):
  fq1=sample_table.loc[ wildcards.sample, "fq1" ]
  fq2=sample_table.loc[ wildcards.sample, "fq2" ]
  return {"r1": fq1, "r2": fq2 }

# define a function for getting the read group information
# from the sample table for each particular sample (according
# to the wildcard value)
def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}_{sample_id}_{library}_{flowcell}_{lane}\tSM:{sample_id}\tPL:{platform}\tLB:{library}\tPU:{flowcell}.{lane}.{library}'".format(
        sample=wildcards.sample,
        sample_id=sample_table.loc[(wildcards.sample), "sample_id"],
        platform=sample_table.loc[(wildcards.sample), "platform"],
        library=sample_table.loc[(wildcards.sample), "library"],
        flowcell=sample_table.loc[(wildcards.sample), "flowcell"],
        lane=sample_table.loc[(wildcards.sample), "lane"],
    )




### Specify rule "all"
# By default, Snakemake tries to create the input files needed
# for the first rule in the Snakefile, so we define the first
# rule to ask for results/vcf/all.vcf.gz
rule all:
  input: 
    "results/vcf-stats/all.vcf.stats"





rule genome_faidx:
  input:
    "data/genome/OmykA.fasta",
  output:
    "data/genome/OmykA.fasta.fai",
  conda:
    "envs/bwa2sam.yaml"
  log:
    "results/logs/genome_faidx.log",
  shell:
    "samtools faidx {input} 2> {log} "


rule genome_dict:
  input:
    "data/genome/OmykA.fasta",
  output:
    "data/genome/OmykA.dict",
  conda:
    "envs/bwa2sam.yaml"
  log:
    "results/logs/genome_dict.log",
  shell:
    "samtools dict {input} > {output} 2> {log} "


rule bwa_index:
  input:
    "data/genome/OmykA.fasta"
  output:
    multiext("data/genome/OmykA.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
  conda:
    "envs/bwa2sam.yaml"
  log:
    out="results/logs/bwa_index/bwa_index.log",
    err="results/logs/bwa_index/bwa_index.err"
  shell:
    "bwa-mem2 index {input} > {log.out} 2> {log.err} "




rule trim_reads:
  input:
    unpack(get_fastqs) # unpack creates named inputs from the dict that
                       # get_fastqs returns
  output:
    r1="results/trimmed/{sample}_R1.fastq.gz",
    r2="results/trimmed/{sample}_R2.fastq.gz",
    html="results/qc/fastp/{sample}.html",
    json="results/qc/fastp/{sample}.json"
  conda:
    "envs/fastp.yaml"
  log:
    out="results/logs/trim_reads/{sample}.log",
    err="results/logs/trim_reads/{sample}.err",
  params:
    as1=config["params"]["fastp"]["adapter_sequence1"],
    as2=config["params"]["fastp"]["adapter_sequence2"],
    parm=config["params"]["fastp"]["other_options"]
  shell:
    " fastp -i {input.r1} -I {input.r2}       "
    "       -o {output.r1} -O {output.r2}     "
    "       -h {output.html} -j {output.json} "
    "  --adapter_sequence={params.as1}        "
    "  --adapter_sequence_r2={params.as2}     "
    "  {params.parm} > {log.out} 2> {log.err} "
    


rule map_reads:
  input:
    r1="results/trimmed/{sample}_R1.fastq.gz",
    r2="results/trimmed/{sample}_R2.fastq.gz",
    genome="data/genome/genome.fasta",
    idx=multiext("data/genome/genome.fasta", ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
  output:
    "results/bam/{sample}.bam"
  conda:
    "envs/bwa2sam.yaml"
  log:
    "results/logs/map_reads/{sample}.log"
  threads: 2
  resources:
    mem_mb=7480,
    time="01:00:00"
  params:
    RG=get_read_group
  shell:
    " (bwa-mem2 mem -t {threads} {params.RG} {input.genome} {input.r1} {input.r2} | "
    " samtools view -u | "
    " samtools sort - > {output}) 2> {log} "




rule mark_duplicates:
  input:
    "results/bam/{sample}.bam"
  output:
    bam="results/mkdup/{sample}.bam",
    bai="results/mkdup/{sample}.bai",
    metrics="results/qc/mkdup_metrics/{sample}.metrics"
  conda:
    "envs/gatk.yaml"
  log:
    "results/logs/mark_duplicates/{sample}.log"
  shell:
    " gatk MarkDuplicates  "
    "  --CREATE_INDEX "
    "  -I {input} "
    "  -O {output.bam} "
    "  -M {output.metrics} > {log} 2>&1 "




rule make_gvcfs_by_chromo:
  input:
    bam="results/mkdup/{sample}.bam",
    bai="results/mkdup/{sample}.bai",
    ref="data/genome/genome.fasta",
    idx="data/genome/genome.dict",
    fai="data/genome/genome.fasta.fai"
  output:
    gvcf="results/gvcf/{chromo}/{sample}.g.vcf.gz",
    idx="results/gvcf/{chromo}/{sample}.g.vcf.gz.tbi",
  conda:
    "envs/gatk.yaml"
  log:
    "results/logs/make_gvcfs_by_chromo/{chromo}/{sample}.log"
  params:
    java_opts="-Xmx4g",
    hmt=config["params"]["gatk"]["HaplotypeCaller"]["hmm_threads"]
  shell:
    " gatk --java-options \"{params.java_opts}\" HaplotypeCaller "
    " -R {input.ref} "
    " -I {input.bam} "
    " -O {output.gvcf} "
    " -L {wildcards.chromo}    "           
    " --native-pair-hmm-threads {params.hmt} " 
    " -ERC GVCF > {log} 2> {log} "




rule import_genomics_db_by_chromo:
  input:
    gvcfs=expand("results/gvcf/{{chromo}}/{s}.g.vcf.gz", s=SAMPLES)
  output:
    gdb=directory("results/genomics_db/{chromo}")
  conda:
    "envs/gatk.yaml"
  log:
    "results/logs/import_genomics_db_by_chromo/{chromo}.log"
  params:
    java_opts="-Xmx4g"
  shell:
    " VS=$(for i in {input.gvcfs}; do echo -V $i; done); "  # make a string like -V file1 -V file2
    " gatk --java-options \"-Xmx4g\" GenomicsDBImport "
    "  $VS  "
    "  --genomicsdb-workspace-path {output.gdb} "
    "  -L  {wildcards.chromo} 2> {log} "




rule vcf_from_gdb_by_chromo:
  input:
    gdb="results/genomics_db/{chromo}",
    ref="data/genome/genome.fasta",
    fai="data/genome/genome.fasta.fai",
    idx="data/genome/genome.dict",
  output:
    vcf="results/chromo_vcfs/{chromo}.vcf.gz",
    idx="results/chromo_vcfs/{chromo}.vcf.gz.tbi",
  conda:
    "envs/gatk.yaml"
  log:
    "results/logs/vcf_from_gdb_by_chromo/{chromo}.txt"
  shell:
    " gatk --java-options \"-Xmx4g\" GenotypeGVCFs "
    "  -R {input.ref}  "
    "  -V gendb://{input.gdb} "
    "  -O {output.vcf} 2> {log} "



# break out just the indels
rule select_indels:
  input:
    vcf="results/chromo_vcfs/{chromo}.vcf.gz",
    idx="results/chromo_vcfs/{chromo}.vcf.gz.tbi"
  output:
    "results/hard_filtering/indels/{chromo}.vcf.gz"
  conda:
    "envs/gatk.yaml"
  log:
    "results/logs/select_indels/{chromo}.log"
  shell:
    "gatk SelectVariants "
    " -V {input.vcf} "
    " -select-type INDEL "
    " -O {output} > {log} 2>&1 " 



# break out just the snps
rule select_snps:
  input:
    vcf="results/chromo_vcfs/{chromo}.vcf.gz",
    idx="results/chromo_vcfs/{chromo}.vcf.gz.tbi"
  output:
    "results/hard_filtering/snps/{chromo}.vcf.gz"
  conda:
    "envs/gatk.yaml"
  log:
    "results/logs/select_indels/{chromo}.log"
  shell:
    "gatk SelectVariants "
    " -V {input.vcf} "
    " -select-type SNP "
    " -O {output} > {log} 2>&1 "


# filter the indels according to the recommendations at:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
rule hard_filter_indels:
  input:
    vcf="results/hard_filtering/indels/{chromo}.vcf.gz"
  output:
    vcf="results/hard_filtering/indels-filtered/{chromo}.vcf.gz"
  conda:
    "envs/gatk.yaml"
  log:
    "results/logs/hard_filter_indels/{chromo}.log"
  params:
    filters=config["params"]["gatk"]["VariantFiltration"]["indels"]
  shell:
    "gatk VariantFiltration "
    " -V {input.vcf} "
    " {params.filters}      "
    " -O {output} > {log} 2>&1 "


# filter the snps according to the recommendations at:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-fi
rule hard_filter_snps:
  input:
    vcf="results/hard_filtering/snps/{chromo}.vcf.gz"
  output:
    vcf="results/hard_filtering/snps-filtered/{chromo}.vcf.gz"
  conda:
    "envs/gatk.yaml"
  log:
    "results/logs/hard_filter_snps/{chromo}.log"
  params:
    filters=config["params"]["gatk"]["VariantFiltration"]["snps"]
  shell:
    "gatk VariantFiltration "
    " -V {input.vcf} "
    " {params.filters}      "
    " -O {output} > {log} 2>&1 "




rule concat_vcfs:
  input:
    vcfs=expand("results/missing-corrected/{c}.vcf.gz", c=CHROMOS)
  output:
    vcf="results/vcf/all.vcf.gz"
  conda:
    "envs/bcftools.yaml"
  log:
    "results/concat_vcfs/all.log"
  shell:
    "bcftools concat -n {input.vcfs} > {output.vcf} 2> {log} "




rule merge_vcfs:
  input:
    indels="results/hard_filtering/indels-filtered/{chromo}.vcf.gz",
    snps="results/hard_filtering/snps-filtered/{chromo}.vcf.gz"
  output:
    vcf="results/hard_filtering/merged/{chromo}.vcf.gz"
  conda:
    "envs/gatk.yaml"
  log:
    "results/merge_vcfs/{chromo}.log"
  benchmark:
    "benchmarks/merge_vcfs/{chromo}.tsv"
  shell:
    "gatk MergeVcfs -I {input.indels} -I {input.snps} -O {output.vcf} 2> {log} "




rule correct_merged_vcfs:
  input:
    vcfs="results/hard_filtering/merged/{chromo}.vcf.gz"
  output:
    vcfs="results/missing-corrected/{chromo}.vcf.gz"
  conda:
    "envs/bcftools.yaml"
  log:
    "results/correct_merged_vcfs/{chromo}.log"
  shell:
    "(bcftools +setGT {input.vcfs} -- -t q -n . -i 'FMT/DP=0 | (FMT/PL[:0]=0 & FMT/PL[:1]=0 & FMT/PL[:2]=0)' | "
    "bcftools +fill-tags - -- -t 'NMISS=N_MISSING' | "
    "bcftools view -Oz - > {output.vcfs}; "
    "bcftools index -t {output.vcfs}) 2> {log} "




rule vcf_stats:
  input:
    "results/vcf/all.vcf.gz"
  output:
    "results/vcf-stats/all.vcf.stats"
  conda:
    "envs/bcftools.yaml"
  log:
    "results/vcf-stats/all.log"
  shell:
    "bcftools stats {input} > {output} 2> {log}"
