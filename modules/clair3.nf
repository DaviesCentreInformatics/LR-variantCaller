/**
 * CLAIR3 process for calling small variants (SNPs and small indels) from long reads
 * Uses deep learning models trained specifically for long-read data
 */
process CLAIR3 {
	tag "$sampleID"
	label "process_high","error_retry"

	publishDir "$params.outdir/variants/SNPs/$sampleID/$ctg_name", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)  // Aligned reads in BAM format with index
	path fa                                    // Reference genome FASTA
	path faidx                                 // Reference genome FASTA index

	output:
	path("*.vcf.gz")         , emit: vcf        // Compressed VCF with variant calls
	path("*.vcf.gz.tbi")     , emit: vcf_tbi    // VCF index
	tuple val(sampleID), path("*.gvcf.gz")        , emit: gvcf       // Genomic VCF for future joint calling
	path("*.gvcf.gz.tbi")    , emit: gvcf_tbi   // GVCF index
	path("log")              , emit: log_dir    // Log directory
	path("*.log")            , emit: log        // Log file
	


	script:
	//model_path = "/opt/models/r941_prom_hac_g238/"
	model_path = params.model_path                // Path to Clair3 model (inside the Clair3 container)
	bam_path = bam.toString()
	//ctg_pattern = "$params.chrom_pattern"
	//ctg_name = bam_path =~ ctg_pattern
	bam_name = bam_path.find(/(?<=_)(\d{1,2}|X|Y)(?=\.bam)/)  // Extract chromosome name from BAM filename
	//ctg_name = ctg_name[0]
	ctg_name = bam_name
	
	// work_dir="${params.temp_dir}"
	// echo "\${work_dir}"
	// clair3_temp=`mktemp -d \${work_dir}/CLAIRXXXXXX`

	"""
	# Run Clair3 variant caller for SNPs and small indels
	run_clair3.sh --bam_fn ${bam} \
		--ref_fn=${fa} \
		--threads=${task.cpus} \
		--platform=ont \
		--sample_name=${sampleID} \
		--model_path=${model_path} \
		--ctg_name=${ctg_name} \
		--gvcf \
		--output=./
	"""
	// mkdir -p $params.outdir/variants/SNPs/$sampleID/$ctg_name
	// cp -rf \${clair3_temp}/* $params.outdir/variants/SNPs/$sampleID/$ctg_name
}