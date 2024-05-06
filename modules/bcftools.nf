process BCFTOOLS_STATS {
	tag "$sampleID"
	label "process_low"

	publishDir "$params.outdir/bcftools/stats", mode: 'copy'

	input:
	tuple val(sampleID), path(vcf)

	output:
	path "*.txt", emit: stats

	script:
	outname = "${vcf.baseName}".replace(".vcf.gz", "")
	"""
	if [[ ! "${vcf}" =~ "_MT" ]]; then
		bcftools stats --threads ${task.cpus} \
		${vcf} > ${outname}.txt
	
	else
		echo "Skipping ${sampleID} MT file." > ${outname}.txt
	fi
	"""
}