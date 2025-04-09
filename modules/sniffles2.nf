/**
 * SNIFFLES2 process for detecting structural variants from long read alignments
 * Sniffles2 is optimized for calling SVs from ONT or PacBio reads
 * with improved performance over the original Sniffles
 */
process SNIFFLES2 {
	tag "$sampleID"
	label "process_medium", "error_retry"

	publishDir "$params.outdir/variants/SVs/sniffles/$sampleID", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)  // Aligned reads in BAM format with index
	path fa                                    // Reference genome FASTA
	path faidx                                 // Reference genome FASTA index

	output:
	path "*.snf.bz2", emit: snf                // Sniffles snf output format (compressed)
	path "*.sniffles.vcf.gz", emit: snf_vcf    // VCF file with SV calls (compressed)
	tuple val(sampleID), path("*.sniffles.vcf.gz"), emit: res_tuple  // Sample with VCF for downstream processes
	// Need to create a more specific output definition

	// TODO: Make this pattern more generalisable or configurable.
	script:
	bam_path = bam.toString()
	bam_name = bam_path.find(/(?<=_)(\d{1,2}|X|Y|MT|M)(?=\.bam)/)  // Extract chromosome name from BAM filename
	//ctg_name = ctg_name[0]
	//ctg_name = bam_name

	"""
	# Run Sniffles2 with parameters
	sniffles --input $bam \
		--vcf ${sampleID}.sniffles.vcf.gz \
		--snf ${sampleID}.snf \
		--reference $fa \
        --minsvlen 50 \
        --mapq 20 \
		--threads ${task.cpus} \
		--output-rnames
	bzip2 ${sampleID}.snf
	"""
}