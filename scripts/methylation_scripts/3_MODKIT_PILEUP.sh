#!/bin/bash -l
#SBATCH -p a100cpu,icelake                                         # partition (this is the queue your job will be added to)
#SBATCH -n 1                                                    # number of tasks (sequential job starts 1 task) (check this if your job unexpectedly uses 2 nodes)
#SBATCH -c 16                                                    # number of cores (sequential job calls a multi-thread program that uses 8 cores)
#SBATCH --time=0-05:00                                          # time allocation, which has the format (D-HH:MM), here set to 1 hour
#SBATCH --mem=32GB                                              # specify memory required per node (here set to 16 GB)

# Configure notifications 
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=callum.macphillamy@adelaide.edu.au          # Email to which notifications will be sent


# This script performs methylation analysis using Modkit and SAMtools.
# It is designed to be executed on a SLURM-based high-performance computing (HPC) cluster.
#
# Usage:
#   ./3_MODKIT_PILEUP.sh <bam> <ref> <outpath>
#
# Arguments:
#   <bam>      : Input BAM file containing aligned sequencing reads.
#   <ref>      : Reference genome file in FASTA format.
#   <outpath>  : Output directory where results will be saved.
#
# SLURM Configuration:
#   - Partition: a100cpu, icelake
#   - Tasks: 1 (sequential job)
#   - Cores: 16 (multi-threaded processing)
#   - Time: 5 hours
#   - Memory: 32GB
#   - Notifications: Email on job completion or failure
#   - Email: callum.macphillamy@adelaide.edu.au
#
# Modules Required:
#   - HTSlib: For handling high-throughput sequencing data.
#   - SAMtools: For indexing and processing BAM files.
#
# Workflow:
#   1. Index the input BAM file using SAMtools with 16 threads.
#   2. Run Modkit's pileup command to generate methylation data:
#      - Logs are saved to a file named after the input BAM file.
#      - Output is a BED file containing methylation information.
#
# Output:
#   - Log file: <outpath>/<basename>.modkit.pileup.log
#   - BED file: <outpath>/<basename>.bedMethyl


module load HTSlib
module load SAMtools

bam=$1
ref=$2
outpath=$3

samtools index -@ 16 ${bam}
/hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/SOFTWARE/modkit-0.3.0/modkit pileup \
        --log-filepath ${outpath}/$(basename -s .bam $bam).modkit.pileup.log \
        -t 16 \
        -r ${ref} \
        --preset traditional ${bam} ${outpath}/$(basename -s .bam $bam).bedMethyl