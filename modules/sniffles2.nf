process SNIFFLES2 {
	tag "$sampleID"
	label "process_medium", "error_retry"

	publishDir "$params.outdir/variants/SVs/$sampleID", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)
	path fa
	path faidx

	output:
	path "*"
	// Need to create a more specific output definition

	// TODO: Make this pattern more generalisable or configurable.
	script:
	bam_path = bam.toString()
	bam_name = bam_path.find(/(?<=_).*(?=\.bam)/)
	//ctg_name = ctg_name[0]
	ctg_name = bam_name

	"""
	sniffles --input $bam \
		--vcf ${sampleID}.${ctg_name}.vcf.gz \
		--snf ${sampleID}.${ctg_name}.snf \
		--reference $fa \
		--threads ${task.cpus} \
		--output-rnames
	bzip2 ${sampleID}.${ctg_name}.snf
	"""
}