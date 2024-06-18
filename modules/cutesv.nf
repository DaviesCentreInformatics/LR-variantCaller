process CUTESV {
	tag "$sampleID"
	label "process_high","error_retry"

	publishDir "$params.outdir/variants/SVs/cuteSV/$sampleID", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)
	path fa 
	path faidx

	output:
	path "*.cuteSV.vcf.gz", emit: cu_vcf
	tuple val(sampleID), path("*.vcf.gz"), emit: res_tuple

	script:
	bam_path = bam.toString()
	bam_name = bam_path.find(/(?<=_)(\d{1,2}|X|Y|MT|M)(?=\.bam)/)
	ctg_name = bam_name
	"""
	cuteSV --threads ${task.cpus} \
		--max_cluster_bias_INS 100 \
		--max_cluster_bias_DEL 100 \
		--diff_ratio_merging_INS 0.3 \
		--diff_ratio_merging_DEL 0.3 \
		--min_support 5 \
		--min_size 50 \
		--genotype \
		-S ${sampleID} \
		--report_readid \
		${bam} \
		${fa} \
		${sampleID}_${ctg_name}.cuteSV.vcf \
		.

	bzip2 ${sampleID}_${ctg_name}.cuteSV.vcf
	"""
}