process DYSGU {
	tag "$sampleID"
	label "process_medium", "error_retry"

	publishDir "$params.outdir/variants/SVs/dysgu/$sampleID", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)
	path fa
	path faidx

	output:
	path "*.vcf.gz", emit: dys_vcf
	tuple val(sampleID), path("*.dysgu.vcf.gz"), emit: res_tuple
	
	script:
	bam_path = bam.toString()
	bam_name = bam_path.find(/(?<=_)(\d{1,2}|X|Y|MT|M)(?=\.bam)/)
	//ctg_name = ctg_name[0]
	ctg_name = bam_name

	"""
	temp_dir="/tmp/\$RANDOM"

	dysgu call --mode nanopore \
		-p ${task.cpus} \
		-c \
		${fa} \
		\$temp_dir \
		${bam} | bzip2 - > ${sampleID}_${ctg_name}.dysgu.vcf.gz
	"""
}