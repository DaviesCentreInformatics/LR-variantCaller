/**
 * FILTLONG process to filter and trim long reads
 * Selects reads based on quality, length, and other metrics
 */
process FILTLONG {
	tag "$sampleID"
	label 'process_low', 'error_retry'
	
	//debug true

	publishDir "$params.outdir/filtlong", mode: 'copy',
		saveAs: { filename -> 
			if (filename.endsWith('.log') || filename.endsWith('.err')) "logs/${filename}"
			else null }

	input:
	tuple val(sampleID), path(bam)  // Sample ID and input reads in BAM format

	output:
	tuple val(sampleID), path("${sampleID}.filtered.bam") , emit: result_tuple  // Filtered reads in BAM format
	// path "${sampleID}.trimmed.fastq.gz"                       , emit: trimmed_fastq
	// path "${sampleID}.filtlong.log"                           , emit: log
	// path "${sampleID}.filtlong.err"                           , emit: err
	// path "version.txt"                                        , emit: version
	
	script:
	"""
	# Convert BAM to FASTQ for filtlong processing
	echo "Converting ${bam} to fastq..." 
	samtools fastq --verbosity 3 -@ ${task.cpus} ${bam} | pigz -p ${task.cpus} - > ${sampleID}.fastq.gz
	
	# Filter reads using filtlong (min length 200)
	echo "Filtering ${sampleID}.fastq..."
	filtlong --min_length 200 ${sampleID}.fastq.gz | \
	seqkit seq -j ${task.cpus} -n | cut -f1 -d ' ' \
		> ${sampleID}.filtered_readnames.txt
	
	# Extract filtered reads from original BAM
	echo "Taking filtered reads from ${bam}..."
	samtools view -@ ${task.cpus} -b \
		-N ${sampleID}.filtered_readnames.txt ${bam} \
		> ${sampleID}.filtered.bam
	"""
}