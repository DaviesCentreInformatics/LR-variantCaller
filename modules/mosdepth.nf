process MOSDEPTH {
	tag "$sampleID"
	label "process_medium", "error_retry"

	publishDir "$params.outdir/mosdepth", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)

	output:
	path "*", emit: report

	script:
	"""
	mosdepth -t ${task.cpus} \
		--by 1000 \
		--no-per-base \
		$sampleID \
		$bam
	"""
}