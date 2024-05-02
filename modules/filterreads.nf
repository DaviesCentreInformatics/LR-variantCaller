process FILTER_READS {
	tag "$sampleID"
	label 'process_low', 'error_retry'
	
	//debug true

	publishDir "$params.outdir/filtered_reads", mode: 'copy',
		saveAs: { filename -> 
			if (filename.endsWith('.log') || filename.endsWith('.err')) "logs/${filename}"
			else filename }

	input:
	tuple val(sampleID), path(fastq)
	tuple val(sampleID), path(bam)

	output:
	tuple val(sampleID), path("${sampleID}.filtered.bam") , emit: result_tuple
	
	script:
	"""
	seqkit seq -j ${task.cpus} -n ${fastq} | cut -f1 -d' ' > ${sampleID}.filtered_readnames.txt
	samtools view -@ ${task.cpus} -b -N ${sampleID}.filtered_readnames.txt ${bam} > ${sampleID}.filtered.bam
	"""
}