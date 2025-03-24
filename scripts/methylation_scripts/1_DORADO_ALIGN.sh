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
