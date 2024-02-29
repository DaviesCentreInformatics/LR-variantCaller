process NANOPLOT {
	tag "$sampleID"
	label "process_low", "error_retry"

	publishDir "$params.outdir/nanoplot", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("trimmed") > 0 ) "trimmed/$filename"
			else "raw/$filename" }

	input:
	tuple val(sampleID), path(fastq)

	output:
	path "${sampleID}${fastq.toString().indexOf("trimmed") > 0 ? "_trimmed" : "_raw"}/*", emit: report

	script:
	outdir = "${sampleID}${fastq.toString().indexOf("trimmed") > 0 ? "_trimmed" : "_raw"}"
	outprefix = "${sampleID}${fastq.toString().indexOf("trimmed") > 0 ? "_trimmed_" : "_raw_"}"

	"""
	NanoPlot -t ${task.cpus} \
		--fastq ${fastq} \
		--outdir $outdir \
		-p $outprefix \
		--loglength \
		--plots dot
	"""
}