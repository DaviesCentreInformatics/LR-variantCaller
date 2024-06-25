#!/usr/bin/env bash

bam=$1
outdir=$2
sampleID=$3

echo "Converting ${bam} to fastq..." \
samtools fastq --verbosity 3 -@ 8 ${bam} | pigz -p 8 - > ${outdir}/${sampleID}.fastq.gz \
	
echo "Filtering ${sampleID}.fastq..." \
filtlong --min_length 200 ${outdir}/${sampleID}.fastq.gz | \
seqkit seq -j 10 -n | cut -f1 -d ' ' \
		> ${outdir}/${sampleID}.filtered_readnames.txt

echo "Taking filtered reads from ${bam}..." \
samtools view -@ 16 -b \
		-N ${outdir}/${sampleID}.filtered_readnames.txt ${bam} \
		> ${outdir}/${sampleID}.filtered.bam