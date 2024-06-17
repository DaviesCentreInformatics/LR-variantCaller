process MODKIT {
    tag "$sampleID"
	label "process_medium", "error_retry"

	debug true

	publishDir "$params.outdir/modkit/$sampleID", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)
	path fa
	path faidx

	output:
	path "*.bedMethyl", emit: bedMethyl
	path "*.modkit.log", emit: modkit_log
	
	script:
	bam_path = bam.toString()
	bam_name = bam_path.find(/(?<=_)(\d{1,2}|X|Y|MT|M)(?=\.bam)/)
	//ctg_name = ctg_name[0]
	ctg_name = bam_name

	"""
	modkit pileup --log-filepath ${sampleID}.${ctg_name}.modkit.log \
		--region ${ctg_name} \
		-t ${task.cpus} \
		--ref ${fa} \
		--preset traditional \
		${bam} \
		${sampleID}.${ctg_name}.bedMethyl
	"""
}