/**
 * CUTESV process for detecting structural variants from long read alignments
 * Detects deletions, insertions, inversions, duplications, and translocations
 */
process CUTESV {
	tag "$sampleID"
	label "process_high","error_retry"

	publishDir "$params.outdir/variants/SVs/cuteSV/$sampleID", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)  // Aligned reads in BAM format with index
	path fa                                    // Reference genome FASTA
	path faidx                                 // Reference genome FASTA index

	output:
	path "*.cuteSV.vcf.gz", emit: cu_vcf                  // Compressed VCF with structural variants
	tuple val(sampleID), path("*.vcf.gz"), emit: res_tuple  // Sample ID with VCF for downstream processes

	script:
	// bam_path = bam.toString()
	// bam_name = bam_path.find(/(?<=_)(\d{1,2}|X|Y|MT|M)(?=\.bam)/)
	// ctg_name = bam_name
	"""
	# Run CuteSV with specific parameters for structural variant detection
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
		${sampleID}.cuteSV.vcf \
		.

	gzip ${sampleID}.cuteSV.vcf
	"""
}