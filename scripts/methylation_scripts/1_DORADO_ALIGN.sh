#!/bin/bash -l

# Configure the resources required
#SBATCH -p a100cpu,icelake 
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --time=0-15:00
#SBATCH --mem=64GB
# Configure notifications 
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=callum.macphillamy@adelaide.edu.au          # Email to which notifications will be sent

# This script performs alignment of sequencing reads to a reference genome using the Dorado aligner.
# It requires three input arguments: the reference genome file, the reads file, and the output file name.
#
# Usage:
#   ./1_DORADO_ALIGN.sh <reference> <reads> <output>
#
# Arguments:
#   <reference>  Path to the reference genome file.
#   <reads>      Path to the sequencing reads file.
#   <output>     Path to the output BAM file.
#
# Modules:
#   - HTSlib: Required for handling high-throughput sequencing data.
#   - SAMtools: Used for converting the alignment output to BAM format.
#
# Notes:
#   - The script uses the Dorado aligner with verbose mode enabled (-v).
#   - It processes the alignment using a single node (-N 1) and 32 threads (-t 32).
#   - The output is piped directly to SAMtools for conversion to BAM format.
#   - Ensure that the Dorado binary path is correct and accessible.
module load HTSlib
module load SAMtools

reference=$1
reads=$2
output=$3

/hpcfs/users/a1767591/DORADO/dorado-0.5.3-linux-x64/bin/dorado aligner -v \
        -N 1 \
        -t 32 \
        $reference \
        $reads | samtools view -b -@ 32 > $output
