/**
 * MULTIQC process to aggregate quality control metrics
 * Combines various QC reports into a single interactive HTML report
 */
process MULTIQC {
	tag "multiqc"
	label "process_medium"


	publishDir "$params.outdir/multiqc", mode: 'copy'

	input:
	path multiqc_files, stageAs: "?/*"  // All QC files to be aggregated

	output:
	path "*multiqc_report.html"  // HTML report
	path "multiqc_data"          // Data directory with parsed statistics

	script:
	"""
	# Run MultiQC on all input files in the current directory
	multiqc .
	"""
}