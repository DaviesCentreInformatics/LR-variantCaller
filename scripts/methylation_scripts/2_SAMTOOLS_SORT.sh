#!/bin/bash -l


# Configure the resources required
#SBATCH -p a100cpu,icelake                                         # partition (this is the queue your job will be added to)
#SBATCH -n 1                                                    # number of tasks (sequential job starts 1 task) (check this if your job unexpectedly uses 2 nodes)
#SBATCH -c 32                                                    # number of cores (sequential job calls a multi-thread program that uses 8 cores)
#SBATCH --time=0-03:00                                          # time allocation, which has the format (D-HH:MM), here set to 1 hour
#SBATCH --mem=32GB                                              # specify memory required per node (here set to 16 GB)

# Configure notifications 
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=callum.macphillamy@adelaide.edu.au          # Email to which notifications will be sent


# This script is a SLURM batch job script for sorting and indexing BAM files using SAMtools.
# It is designed to be executed on a high-performance computing cluster with SLURM workload manager.

# SLURM Directives:
# - Specifies the partition(s) to use for the job (`a100cpu,icelake`).
# - Allocates 1 task (`-n 1`) and 32 CPU cores (`-c 32`) for multi-threaded processing.
# - Sets a time limit of 3 hours (`--time=0-03:00`) and allocates 32GB of memory (`--mem=32GB`).
# - Configures email notifications to be sent on job completion or failure to the specified email address.

# Modules:
# - Loads the required modules `HTSlib` and `SAMtools` for processing BAM files.

# Input:
# - The script expects a single argument (`$1`), which is the path to the input BAM file.

# Processing:
# - Uses `samtools sort` to sort the input BAM file with 32 threads and outputs it as `sorted_<input_bam>`.
# - Uses `samtools index` to create an index for the sorted BAM file with 32 threads.

# Output:
# - A sorted BAM file named `sorted_<input_bam>`.
# - An index file for the sorted BAM file.

# Usage:
# - Submit the script to the SLURM scheduler with the input BAM file as an argument:
#   `sbatch 2_SAMTOOLS_SORT.sh <input_bam>`


module load HTSlib
module load SAMtools

bam=$1

samtools sort -@ 32 ${bam} > sorted_${bam}
samtools index -@ 32 sorted_${bam}