/**
 * NANOPLOT process for quality control of long sequencing reads
 * Creates various plots and statistics to assess read quality,
 * read length distribution, and other metrics for ONT data
 */
process NANOPLOT {
	tag "$sampleID"
	label "process_low", "error_retry"

	publishDir "$params.outdir/nanoplot", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("filtered") > 0 ) "filtered/$filename"
			else "raw/$filename" }

	input:
	tuple val(sampleID), path(reads)  // Sample ID and sequencing reads (FASTQ or BAM)

	output:
	path "${sampleID}${reads.toString().indexOf("filtered") > 0 ? "_filtered" : "_raw"}/*", emit: report  // NanoPlot reports and figures

	script:
	// Determine output directory and prefix based on whether reads are filtered
	outdir = "${sampleID}${reads.toString().indexOf("filtered") > 0 ? "_filtered" : "_raw"}"
	outprefix = "${sampleID}${reads.toString().indexOf("filtered") > 0 ? "_filtered_" : "_raw_"}"

	// Determine input file type (FASTQ or BAM)
	hts_type = reads.toString().indexOf("fastq") > 0 ? "--fastq" : "--ubam"

	"""
	# Run NanoPlot on the sequencing reads
	NanoPlot -t ${task.cpus} \
		${hts_type} ${reads} \
		--outdir $outdir \
		-p $outprefix \
		--loglength \
		--plots dot \
		--verbose
	"""
}