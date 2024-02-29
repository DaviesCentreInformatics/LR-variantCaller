process CLAIR3 {
	tag "$sampleID"
	label "process_high","error_retry"

	publishDir "$params.outdir/variants/SNPs/$sampleID/$ctg_name", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)
	path fa 
	path faidx

	output:
	//path("*.vcf.gz")
	//path("*.vcf.gz.tbi")
	//path("*.gvcf.gz")
	//path("*.gvcf.gz.tbi")
	//path("log")
	//path("*.log")
	path "*"


	script:
	model_path = "/opt/models/r941_prom_hac_g238/"
	bam_path = bam.toString()
	ctg_pattern = ~/NC_\d+\.\d+/
	ctg_name = bam_path =~ ctg_pattern
	ctg_name = ctg_name[0]
	"""
	run_clair3.sh --bam_fn ${bam} \
		--ref_fn=${fa} \
		--threads=${task.cpus} \
		--platform=ont \
		--sample_name=${sampleID} \
		--model_path=${model_path} \
		--ctg_name=${ctg_name} \
		--remove_intermediate_dir \
		--gvcf \
		--output=./
	"""
}