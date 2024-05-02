process BCFTOOLS_STATS {
	tag "$sampleID"
	label "process_low"

	publishDir "$params.outdir/bcftools/stats", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("clair3") > 0 ) "clair3/$filename"
			else if (filename.indexOf("sniffles2") > 0 ) "sniffles2/$filename"
			else if (filename.indexOf("svim") > 0 ) "svim/$filename"
			else if (filename.indexOf("cutesv") > 0 ) "cutesv/$filename"
			else if (filename.indexOf("dysgu") > 0) "dysgu/$filename"
			else null }

	input:
	path(vcf)

	output:
	path "${sampleID}.stats", emit: stats


	script:
	vcf_path = 
	"""
	bcftools stats -@ ${task.cpus} ${bam} > ${sampleID}.stats
	"""
}