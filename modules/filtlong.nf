process FILTLONG {
	tag "$sampleID"
	label 'process_low', 'error_retry'
	
	//debug true

	publishDir "$params.outdir/filtlong", mode: 'copy',
		saveAs: { filename -> 
			if (filename.endsWith('.log') || filename.endsWith('.err')) "logs/${filename}"
			else filename }

	input:
	tuple val(sampleID), path(fastq)

	output:
	tuple val(sampleID), path("${sampleID}.trimmed.fastq.gz") , emit: result_tuple
	path "${sampleID}.trimmed.fastq.gz"                       , emit: trimmed_fastq
	path "${sampleID}.filtlong.log"                           , emit: log
	path "${sampleID}.filtlong.err"                           , emit: err
	path "version.txt"                                        , emit: version
	
	script:
	"""
	filtlong --min_length 200 ${fastq} | bgzip > ${sampleID}.trimmed.fastq.gz
	cp .command.log ${sampleID}.filtlong.log
	cp .command.err ${sampleID}.filtlong.err
	filtlong --version > version.txt
	"""
}