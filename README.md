# Davies Informatics SV Calling Pipeline

## Table of Contents

- [Davies Informatics SV Calling Pipeline](#davies-informatics-sv-calling-pipeline)
	- [Table of Contents](#table-of-contents)
	- [Introduction](#introduction)
	- [Requirements](#requirements)
	- [Installation](#installation)
	- [Usage](#usage)
	- [Road Map](#road-map)
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

Ensure you have Nextflow installed. The easiest way to do this is to create a
conda environment and install Nextflow into it.

``` bash
conda create -n nextflow python=3
```

[Back to top](#)

## Usage

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

nextflow run 
```

Once it's running, you can detach from the screen session by pressing `Ctrl + A`
then `D`. You can reattach to the session by running `screen -r nextflow`.


[Back to top](#)

## Road Map

- [ ] Add support for other long read technologies.
- [ ] Add multiple SV callers
- [ ] Add tool to identify consensus SVs.

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