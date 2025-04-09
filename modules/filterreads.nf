/**
 * FILTER_READS process to extract reads from a BAM file based on a FASTQ file
 * Identifies reads in a FASTQ file and extracts them from a corresponding BAM
 */
process FILTER_READS {
	tag "$sampleID"
	label 'process_low', 'error_retry'
	
	//debug true

	publishDir "$params.outdir/filtered_reads", mode: 'copy',
		saveAs: { filename -> 
			if (filename.endsWith('.log') || filename.endsWith('.err')) "logs/${filename}"
			else filename }

	input:
	tuple val(sampleID), path(fastq)  // Sample ID and FASTQ file containing reads to extract
	tuple val(sampleID), path(bam)    // BAM file from which to extract reads

	output:
	tuple val(sampleID), path("${sampleID}.filtered.bam") , emit: result_tuple  // BAM with only extracted reads
	
	script:
	"""
	# Extract read names from FASTQ file
	seqkit seq -j ${task.cpus} -n ${fastq} | cut -f1 -d' ' > ${sampleID}.filtered_readnames.txt
	
	# Extract reads with those names from the BAM file
	samtools view -@ ${task.cpus} -b -N ${sampleID}.filtered_readnames.txt ${bam} > ${sampleID}.filtered.bam
	"""
}