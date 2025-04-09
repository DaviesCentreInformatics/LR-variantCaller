/**
 * DYSGU process for calling structural variants from long reads
 * Specialized SV caller for nanopore data
 */
process DYSGU {
	tag "$sampleID"
	label "process_medium", "error_retry"

	publishDir "$params.outdir/variants/SVs/dysgu/$sampleID", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)  // Aligned reads in BAM format with index
	path fa                                    // Reference genome FASTA
	path faidx                                 // Reference genome FASTA index

	output:
	path "*.vcf.gz", emit: dys_vcf                          // Compressed VCF with variant calls
	tuple val(sampleID), path("*.dysgu.vcf.gz"), emit: res_tuple  // Sample ID with VCF for downstream processes
	
	script:
	// bam_path = bam.toString()
	// bam_name = bam_path.find(/(?<=_)(\d{1,2}|X|Y|MT|M)(?=\.bam)/)
	// //ctg_name = ctg_name[0]
	// ctg_name = bam_name

	"""
	# Create temporary directory for Dysgu
	temp_dir=$params.temp_dir/${sampleID}_dysgu

	# Run Dysgu SV caller in nanopore mode
	dysgu call --mode nanopore \
		-p ${task.cpus} \
		-c \
        -x \
		${fa} \
		\$temp_dir \
		${bam} | gzip - > ${sampleID}.dysgu.vcf.gz
	"""
}