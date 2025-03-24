#!/bin/bash -l
#SBATCH -p a100cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=10:00:00
#SBATCH --mem=32GB

# Notification configuration
#SBATCH --mail-type=ALL
#SBATCH --mail-user=callum.macphillamy@adelaide.edu.au

module load Singularity

conda activate singularity-env

INPUT_DIR=$1
bam=$2
OUTPUT_DIR=$3
sample_name=$(basename -s .sorted.bam ${bam})
ONT_MODEL_PATH='/hpcfs/users/a1767591/CONTAINERS/clair3_models/r1041_e82_400bps_sup_g615'
#CONTIGS_LIST="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,X,Y,MT"
CONTIGS_LIST="29"

mkdir -p ${OUTPUT_DIR}/${sample_name}
SAMPLE_OUT=${OUTPUT_DIR}/${sample_name}

singularity exec \
	--pid \
	-B ${INPUT_DIR},${SAMPLE_OUT},/hpcfs/users/a1767591/CONTAINERS/clair3_models/r1041_e82_400bps_sup_g615 \
	--bind ${TMPDIR}:/tmp --bind ${TMPDIR} \
	/hpcfs/users/a1767591/CONTAINERS/clair3.sif \
	/opt/bin/run_clair3.sh \
	--bam_fn=${INPUT_DIR}/${bam} \
	--ref_fn=${INPUT_DIR}/ARS_UCD_v2.0.fa \
	--threads=16 \
	--model_path=${ONT_MODEL_PATH} \
	--platform="ont" \
	--ctg_name=${CONTIGS_LIST} \
	--gvcf \
	--remove_intermediate_dir \
	--output=${SAMPLE_OUT}