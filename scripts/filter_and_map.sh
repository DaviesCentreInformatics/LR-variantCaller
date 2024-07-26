#!/usr/bin/env bash

set -eou pipefail

bam=$1
reference=$2
sampleid=$(basename -s .sup.bam ${bam})

# filtlong
samtools fastq --verbosity 3 -@ 8 ${bam} | pigz -p 8 - > ${sampleid}.fastq.gz

filtlong --min_length 200 ${sampleid}.fastq.gz | \
	seqkit seq -j 8 -n | cut -f1 -d ' ' > ${sampleid}.ids

samtools view -@ 16 -b -N ${sampleid}.ids ${bam} > ${sampleid}.filtered.bam

# minimap2
dorado aligner -t 8 ${reference} ${sampleid}.filtered.bam | \
	samtools view -@ 8 -b - > ${sampleid}.aligned.bam

# samtools
samtools sort -@ 16 -o ${sampleid}.wagyu.sorted.bam ${sampleid}.aligned.bam

# Check that the previous step completed successfully
if [ ! -f ${sampleid}.wagyu.sorted.bam ]; then
	echo "Error: ${sampleid}.wagyu.sorted.bam not found"
	exit 1
fi
rm ${sampleid}.fastq.gz ${sampleid}.filtered.bam ${sampleid}.aligned.bam
