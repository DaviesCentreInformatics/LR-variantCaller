# General Steps for CENP-A Analysis

1. Trim and map reads to the reference genome.
``` bash
module use /apps/modules/all/

module load fastp/0.23.4-GCC-11.2.0 
module load Bowtie2/2.5.1-GCC-11.2.0



set -eou pipefail


# If not enough arguments are provided, print usage and exit

if [ $# -ne 5 ]; then
	echo "Usage: $0 <Read1> <Read2> <sample_name> <ref> <out_dir>"
	exit 1
fi

Read1=$1
Read2=$2
sample_name=$3
ref=$4
outdir=$5

read_1_basename=$(basename -s .fastq.gz ${Read1})
read_2_basename=$(basename -s .fastq.gz ${Read2})


mkdir -p $outdir

fastp -i $Read1 -I $Read2 \
	-o ${outdir}/${read_1_basename}.fastq.gz \
	-O ${outdir}/${read_2_basename}.fastq.gz \
	--html ${outdir}/${sample_name}_fastp.html \
	--json ${outdir}/${sample_name}_fastp.json

bowtie2 -k 100 \
	-p 16 \
	-x ${ref} \
	-1 ${outdir}/${read_1_basename}.fastq.gz \
	-2 ${outdir}/${read_2_basename}.fastq.gz |\
	samtools view -@ 16 -b - > ${outdir}/${sample_name}.bam

rm ${outdir}/${read_1_basename}.fastq.gz
rm ${outdir}/${read_2_basename}.fastq.gz

echo "${sample_name} mapping complete."
```

2. Generate a BedGraph file from the BAM file
``` bash
module load BEDTools

set -eou pipefail

if [ $# -ne 4 ]; then
        echo "Usage: $0 <bam> <sample_name> <chrom_sizes> <outdir>"
        exit 1
fi

bam=$1
sample_name=$2
chrom_sizes=$3
outdir=$4

cd ${outdir}
#echo "Sorting ${bam} by coordinates"
#samtools sort --verbosity 4 -@ 16 ${bam} > $(basename -s .bam ${bam}).sorted.bam
#samtools index $(basename -s .bam ${bam}).sorted.bam

#samtools index -@ 16 $(basename ${bam})

echo "Splitting ${bam} by chromosome"
touch ${sample_name}.bedpe.bed

for chrom in $(seq 1 29) X; do
        samtools view -@ 8 -bh ${bam} ${chrom} | samtools sort -n -@ 8 -T chrom_${chrom} > ${outdir}/$(basename -s .bam ${bam}).chrom${chrom}.nameSorted.bam
        
        bedtools bamtobed -bedpe -i ${outdir}/$(basename -s .bam ${bam}).chrom${chrom}.nameSorted.bam >> ${outdir}/${sample_name}.bedpe.bed
        rm ${outdir}/$(basename -s .bam ${bam}).chrom${chrom}.nameSorted.bam
        echo "${sample_name} ${chrom} completed successfully."
done

echo "Filtering ${bam} by fragment length"
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${outdir}/${sample_name}.bedpe.bed > ${outdir}/${sample_name}.clean.bed
cut -f 1,2,6 ${outdir}/${sample_name}.clean.bed | sort -T /hpcfs/users/a1767591/CENPA_analysis/sort_temp -k1,1 -k2,2n -k3,3n > ${outdir}/${sample_name}.fragments.bed
echo "Making final bedGraph for ${sample_name}"
bedtools genomecov -bg -i ${outdir}/${sample_name}.fragments.bed -g ${chrom_sizes} > ${outdir}/${sample_name}.fragments.bedGraph
echo "${sample_name} completed successfully."
```

3. Call peaks with SEACR

``` bash
bash SEACR.sh ${sample_name}.fragments.bedGraph \
    0.1 \
    non \
    stringent \
    ${output}.peaks
```