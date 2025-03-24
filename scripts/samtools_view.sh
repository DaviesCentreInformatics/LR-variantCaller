#!/bin/bash -l
#SBATCH -p a100cpu,skylake,icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=0:00:20
#SBATCH --mem=16GB

# Notification configuration
#SBATCH --mail-type=ALL
#SBATCH --mail-user=callum.macphillamy@adelaide.edu.au

module load SAMtools

bam=$1
chrom=$2
new_bam=$(basename -s .bam ${bam})_${chrom}.bam

samtools view -b ${bam} ${chrom} > ${new_bam}
