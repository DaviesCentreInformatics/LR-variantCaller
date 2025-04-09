/**
 * SVIM process for detecting structural variants from long read alignments
 * This process uses SVIM to identify SVs (insertions, deletions, duplications, inversions)
 * from aligned long reads (e.g., ONT or PacBio)
 */
process SVIM {
	tag "$sampleID"
	label "process_high","error_retry"

	publishDir "$params.outdir/variants/SVs/SVIM/$sampleID", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)  // Aligned reads in BAM format with index
	path fa                                    // Reference genome FASTA
	path faidx                                 // Reference genome FASTA index

	output:
	path "*.svim.vcf.gz", emit: svim_vcf       // Compressed VCF with structural variants
	tuple val(sampleID), path("*.vcf.gz"), emit: res_tuple  // Sample ID with VCF for downstream processes

	script:
	// bam_path = bam.toString()
	// bam_name = bam_path.find(/(?<=_)(\d{1,2}|X|Y|MT|M)(?=\.bam)/)
	// ctg_name = bam_name
	"""
	# Run SVIM with alignment mode on the BAM file
	svim alignment ./ ${bam} ${fa} \
	 	--min_mapq 20 \
		--min_sv_size 50 \
		--minimum_depth 5 \
		--sample ${sampleID} \
		--insertion_sequences
		
	mv variants.vcf ${sampleID}.svim.vcf
	gzip ${sampleID}.svim.vcf
	"""
}