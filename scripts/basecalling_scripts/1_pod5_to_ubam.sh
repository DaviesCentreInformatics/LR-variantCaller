#!/bin/bash -l

# Configure the resources required
#SBATCH -p a100                                                 # partition (this is the queue your job will be added to)
#SBATCH -n 1                                                    # number of tasks (sequential job starts 1 task) (check this if your job unexpectedly uses 2 nodes)
#SBATCH -c 8                                                    # number of cores (sequential job calls a multi-thread program that uses 8 cores)
#SBATCH --time=0-02:00                                          # time allocation, which has the format (D-HH:MM), here set to 1 hour
#SBATCH --gres=gpu:1                                            # generic resource required (here requires 4 GPUs)
#SBATCH --mem=16GB                                              # specify memory required per node (here set to 16 GB)

# Configure notifications 
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=callum.macphillamy@adelaide.edu.au          # Email to which notifications will be sent

module load HTSlib

pod5=$1
out_dir=$2
model=$3
bam_out=$(basename -s .pod5 ${pod5})_${model}.bam

./dorado-0.5.3-linux-x64/bin/dorado basecaller -v \
    -b 0 \
	--min-qscore 9 \
    ./dorado-0.5.3-linux-x64/models/dna_r10.4.1_e8.2_400bps_${model}@v4.2.0 \
    --modified-bases-models ./dorado-0.5.3-linux-x64/models/dna_r10.4.1_e8.2_400bps_${model}@v4.2.0_5mCG_5hmCG@v2/ \
    $pod5 > ${out_dir}/${bam_out}