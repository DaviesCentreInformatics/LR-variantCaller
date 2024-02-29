process MULTIQC {
	tag "multiqc"
	label "process_medium"


	publishDir "$params.outdir/multiqc", mode: 'copy'

	input:
	path multiqc_files, stageAs: "?/*"

	output:
	path "*multiqc_report.html"
	path "multiqc_data"

	script:
	"""
	multiqc .
	"""
}