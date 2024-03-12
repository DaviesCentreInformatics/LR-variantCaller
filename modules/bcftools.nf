process BCFTOOLS_STATS {
	tag "$sampleID"
	label "process_low"

	publishDir "$params.outdir/bcftools", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("flagstats") > 0 ) "flagstats/$filename"
			else null }

	input:
	tuple val(sampleID), path(vcf)
	path fa
	path faidx

	output:
	path "${sampleID}.*.*.stats", emit: vcfstats


	script:
	vcf_path = vcf.toString()
	stats_name = basename(vcf_path).replaceAll(".vcf.gz", "")
	if (vcf_path.indexOf("SV") > 0 ) {

		"""
		bcftools stats -@ ${task.cpus} -F ${fa} ${vcf} > ${stats_name}.SV.stats
		"""
	} else {
		"""
		bcftools stats -@ ${task.cpus} -F ${fa} ${vcf} > ${stats_name}.SNP.stats
		"""
}