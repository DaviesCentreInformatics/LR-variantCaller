process NANOPLOT {
	tag "$sampleID"
	label "process_low", "error_retry"

	publishDir "$params.outdir/nanoplot", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("filtered") > 0 ) "filtered/$filename"
			else "raw/$filename" }

	input:
	tuple val(sampleID), path(reads)

	output:
	path "${sampleID}${reads.toString().indexOf("filtered") > 0 ? "_filtered" : "_raw"}/*", emit: report

	script:
	outdir = "${sampleID}${reads.toString().indexOf("filtered") > 0 ? "_filtered" : "_raw"}"
	outprefix = "${sampleID}${reads.toString().indexOf("filtered") > 0 ? "_filtered_" : "_raw_"}"

	hts_type = reads.toString().indexOf("fastq") > 0 ? "--fastq" : "--ubam"

	"""
	NanoPlot -t ${task.cpus} \
		${hts_type} ${reads} \
		--outdir $outdir \
		-p $outprefix \
		--loglength \
		--plots dot \
		--verbose
	"""
}