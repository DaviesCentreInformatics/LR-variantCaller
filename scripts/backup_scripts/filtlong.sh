#!/bin/bash -l
#SBATCH -p a100cpu,icelake                                      # partition (this is the queue your job will be added to)
#SBATCH -n 1                                                    # number of tasks (sequential job starts 1 task) (check this if your job unexpectedly uses 2 nodes)
#SBATCH -c 16                                                   # number of cores (sequential job calls a multi-thread program that uses 8 cores)
#SBATCH --time=0-10:00                                          # time allocation, which has the format (D-HH:MM), here set to 1 hour
#SBATCH --mem=32GB                                              # specify memory required per node (here set to 16 GB)

# Configure notifications 
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=callum.macphillamy@adelaide.edu.au          # Email to which notifications will be sent

module load Singularity

bam=$1
outdir=$2
sampleID=$(basename -s .sup.bam ${bam})

container="/hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/CONTAINERS/filtlong.sif"

singularity exec --bind ${TMPDIR}:/tmp --bind ${TMPDIR} --bind ${outdir}

echo "Converting ${bam} to fastq..." 
samtools fastq --verbosity 3 -@ 8 ${bam} | pigz -p 8 - > ${outdir}/${sampleID}.fastq.gz
	
echo "Filtering ${sampleID}.fastq..."
filtlong --min_length 200 ${outdir}/${sampleID}.fastq.gz | \
seqkit seq -j 10 -n | cut -f1 -d ' ' \
		> ${outdir}/${sampleID}.filtered_readnames.txt
	
echo "Taking filtered reads from ${bam}..."
samtools view -@ 16 -b \
		-N ${outdir}/${sampleID}.filtered_readnames.txt ${bam} \
		> ${outdir}/${sampleID}.filtered.bam