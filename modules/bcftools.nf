process BCFTOOLS_STATS {
	tag "$sampleID"
	label "process_low"
	debug true

	publishDir "$params.outdir/bcftools", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("SV") > 0 ) "SV/$filename"
			else if (filename.indexOf("SNP") > 0) "SNP/$filename" 
			else null }

	input:
	tuple val(sampleID), path(vcf)
	path fa
	path faidx

	output:
	path "${sampleID}.*.stats", emit: vcfstats


	script:
	vcf_path = vcf.toString().replace(".vcf.gz","")
	stats_name = vcf_path
	if (vcf_path.indexOf("SV") > 0 ) {

		"""
		bcftools stats -@ ${task.cpus} -F ${fa} ${vcf} > ${stats_name}.SV.stats
		"""
	} else {
		"""
		bcftools stats -@ ${task.cpus} -F ${fa} ${vcf} > ${stats_name}.SNP.stats
		"""
	}
}