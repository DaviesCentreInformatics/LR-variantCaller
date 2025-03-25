# Davies Informatics SV Calling Pipeline

## Table of Contents

- [Davies Informatics SV Calling Pipeline](#davies-informatics-sv-calling-pipeline)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Example - Currently being retested with new Sniffles version.](#example---currently-being-retested-with-new-sniffles-version)
  - [Contributing](#contributing)
  - [License](#license)
  - [Acknowledgements](#acknowledgements)
  - [Contact](#contact)

## Introduction

This pipeline was written to call structural variants (SVs) from Oxford Nanopore
Technologies (ONT) long read sequencing data. It is designed to be run on the 
University of Adelaide's HPC, Phoenix. It is written in Nextflow and uses
Singularity to manage the containers.

[Back to top](#)

## Requirements

- [Nextflow](https://www.nextflow.io/)
- [Singularity](https://sylabs.io/guides/3.7/user-guide/installation.html)
- [Docker](https://docs.docker.com/get-docker/)

[Back to top](#)

## Installation

Clone the repository or find it at on Phoenix

``` bash
# If getting from GitHub
git clone git@github.com:DaviesCentreInformatics/LR-variantCaller.git

# If using what's already on Phoenix
ls /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/WORKFLOWS/LR-variantCaller
```

Ensure you have Nextflow installed. The easiest way to do this is to create a
conda environment and install Nextflow into it.

``` bash
conda create -n nextflow python=3
```

[Back to top](#)

## Usage

If you know what you're doing or just need a refresher, you can run the pipeline
with the following command.

``` bash

```bash
nextflow run /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/WORKFLOWS/LR-variantCaller \
    -params-file params.yaml \
    -profile singularity,slurm

# OR
nextflow run /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/WORKFLOWS/LR-variantCaller \
    -params-file params.yaml \
    -profile singularity,slurm \
    -resume
```

## Example - Currently being retested with new Sniffles version.

If this is your first time running the pipeline, this worked example should help
you get started. For the most part, you should be able to copy and paste the
commands as they're written, update the paths and have them work. This example 
uses a Wagyu sample that has been downsampled to ~1x coverage.

> **NOTE:**  
> This fastq file has been downsampled to ~1x coverage. We sequenced to ~10x
> and then `samtools view --subsample 0.1 --subsample-seed 12 sample01.bam > sample01.sub.bam`
> was used to get this file.  
> 10% of ~10x is ~1x.

0. If you haven't already, ensure the reference genome is indexed with 
   `minimap2`. Fortunately, this is very fast compared to `bwa`.

``` bash
minimap2 -d reference.mmi reference.fa
```

We have an already indexed reference available at:  
`/hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/REFERENCES/ARS-UCD_BLRC/ARS_UCD_v2.0.mmi`

1. Create a directory where you want the samplesheet, params file and results to
    be saved.

```bash
# Update this path to reflect your own directory structure.
mkdir /hpcfs/groups/phoenix-hpc-avsci/Callum_MacPhillamy/demo_ont
cd /hpcfs/groups/phoenix-hpc-avsci/Callum_MacPhillamy/demo_ont
mkdir ont_results
mkdir dysgu_temp
```

1. Create a samplesheet. This is a CSV file that contains the sampleID and the 
   full path to the unmapped BAM file. Don't worry that the column header is 
   `fastq`, it's just a placeholder.
   
``` csv
sampleID,fastq
sample01,/hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/TEST_DATA/ONT/sample01.sup.sub.bam
```

This may be useful to make the samplesheet.

```bash
fastq_dir=/hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/TEST_DATA/ONT

(echo "sampleID,fastq" && for f in ${fastq_dir}/*.sup.bam; do echo $( cut -d. -f 1 <(basename $f))","$f; done) > samplesheet.csv

cat samplesheet.csv
# sampleID,fastq
# sample01,/hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/TEST_DATA/ONT/sample01.sup.sub.bam
```

2. Create a params file. This is a YAML file that contains the parameters for 
   the pipeline. You can copy and paste the following into a file called 
   `params.yaml`.

``` yaml
# Clair3 Parameters. These are the paths within the Clair3 container.
model_path: /opt/models/r1041_e82_400bps_sup_v420/

sourceDir: /hpcfs/groups/phoenix-hpc-avsci/Callum_MacPhillamy/demo_ont/ont_results
reference: /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/REFERENCES/ARS-UCD_BLRC/ARS_UCD_v2.0.fa
reference_idx: /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/REFERENCES/ARS-UCD_BLRC/ARS_UCD_v2.0.fa.fai
minimap_index: /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/REFERENCES/ARS-UCD_BLRC/ARS_UCD_v2.0.mmi

# If using mapped bams as input set to true.
is_mapped: false

# If wanting to call SNPs set to true.
call_snps: true

# Set to true, all the SV callers you want to use.
sniffles: true
svim: true
cutesv: true
dysgu: true

temp_dir: /hpcfs/groups/phoenix-hpc-avsci/Callum_MacPhillamy/demo_ont/dysgu_temp
```

Once you have the samplesheet and params file, you can run the pipeline.

``` bash
screen -S sv_caller

module load Singularity

conda activate nextflow

nextflow run /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/WORKFLOWS/LR-variantCaller \
    --samplesheet samplesheet.csv \
    -params-file params.yaml \
    -profile singularity,slurm 
```

Alternatively, you can provide a samplesheet of mapped bam files and call 
variants starting from there.

> **NOTE:**
> You need to ensure the mapped BAM files are sorted and indexed and that these 
> indexes (indices?) are in the same directory as the BAM files.


``` bash
nextflow run /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/WORKFLOWS/LR-variantCaller \
    --samplesheet samplesheet.csv \
	--is_mapped \
	-params-file params.yaml -profile singularity,slurm 
```

Once it's running, you can detach from the screen session by pressing `Ctrl + A`
then `D`. You can reattach to the session by running `screen -r nextflow`.


[Back to top](#)

## Contributing

[Back to top](#)

## License

[GNU-GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

[Back to top](#)

## Acknowledgements

[Back to top](#)

## Contact

Maintainer: Callum MacPhillamy (callum.macphillamy@adelaide.edu.au)

[Back to top](#)