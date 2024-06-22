process SVIM {
	tag "$sampleID"
	label "process_high","error_retry"

	publishDir "$params.outdir/variants/SVs/SVIM/$sampleID", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)
	path fa 
	path faidx

	output:
	path "*.svim.vcf.gz", emit: svim_vcf
	tuple val(sampleID), path("*.vcf.bz2"), emit: res_tuple

	script:
	// bam_path = bam.toString()
	// bam_name = bam_path.find(/(?<=_)(\d{1,2}|X|Y|MT|M)(?=\.bam)/)
	// ctg_name = bam_name
	"""
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