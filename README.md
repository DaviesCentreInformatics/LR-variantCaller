# Davies Informatics ONT SV Calling Pipeline

## Introduction

This pipeline was written to call structural variants (SVs) from Oxford Nanopore
Technologies (ONT) long read sequencing data. It is designed to be run on the 
University of Adelaide's HPC, Phoenix. It is written in Nextflow and uses
Singularity to manage the containers.

## Installation

Ensure you have Nextflow installed. The easiest way to do this is to create a
conda environment and install Nextflow into it.

``` bash
conda create -n nextflow python=3
```

## Setting up for a workflow run

To run the pipeline, first you need to make a samplesheet with the samples you
want to run. The samplesheet should be a CSV file with the following columns:
`sampleID,fastq`. The `sampleID` is the name of the sample and `fastq` is the
complete path to the fastq file for that sample. For example:

``` csv
sampleID,fastq
SRR12898293,/Users/callummacphillamy/Projects/LR-variantCaller/input_data/raw_fastq/SRR12898293.fastq.gz
SRR12898316,/Users/callummacphillamy/Projects/LR-variantCaller/input_data/raw_fastq/SRR12898316.fastq.gz
SRR24678051,/Users/callummacphillamy/Projects/LR-variantCaller/input_data/raw_fastq/SRR24678051.fastq.gz
```

You should create this file in a directory where you want your results to be
saved. This is a requirement for the workflow to work with `singularity`. 
E.g.

``` bash
mkdir -p /hpcfs/groups/phoenix-hpc-avsci/Callum_MacPhillamy/Wagyu_SV/results_batch1
cd /hpcfs/groups/phoenix-hpc-avsci/Callum_MacPhillamy/Wagyu_SV
vim samplesheet_batch1.csv
```

Once you have your samplesheet, we can run the pipeline.

``` bash
screen -S nextflow

module load Singularity

nextflow run /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/WORKFLOWS/LR-variantCaller --samplesheet human_ont_5X_samplesheet.csv --reference $PWD/reference/GCF_009914755.1_T2T-CHM13v2.0_genomic.fa --minimap_index $PWD/reference/GCF_009914755.1_T2T-CHM13v2.0_genomic.mmi --sourceDir /hpcfs/groups/phoenix-hpc-avsci/Callum_MacPhillamy/Wagyu_SV/results_batch1 --outdir /hpcfs/groups/phoenix-hpc-avsci/Callum_MacPhillamy/Wagyu_SV/results_batch1 -profile singularity,slurm
```

Once it's running, you can detach from the screen session by pressing `Ctrl + A`
then `D`. You can reattach to the session by running `screen -r nextflow`.
