#!/bin/bash -l

# Configure the resources required
#SBATCH -p a100cpu,icelake                                      # partition (this is the queue your job will be added to)
#SBATCH -n 1                                                    # number of tasks (sequential job starts 1 task) (check this if your job unexpectedly uses 2 nodes)
#SBATCH -c 16                                                    # number of cores (sequential job calls a multi-thread program that uses 8 cores)
#SBATCH --time=0-04:00                                          # time allocation, which has the format (D-HH:MM), here set to 1 hour
#SBATCH --mem=32GB                                              # specify memory required per node (here set to 16 GB)

# Configure notifications 
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=callum.macphillamy@adelaide.edu.au          # Email to which notifications will be sent

module load HTSlib

bam=$1
out_dir=$2

./dorado-0.6.0-linux-x64/bin/dorado demux -t 16 \
    -v \
    --kit-name "SQK-NBD114-24" \
    --output-dir ${out_dir} \
    ${bam}