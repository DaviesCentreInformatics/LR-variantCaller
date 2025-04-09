/**
 * BCFTOOLS_STATS process to generate statistics from VCF files
 * Provides summary information about variants in a VCF file
 */
process BCFTOOLS_STATS {
	tag "$sampleID"
	label "process_low"

	publishDir "$params.outdir/bcftools/stats", mode: 'copy'

	input:
	tuple val(sampleID), path(vcf)  // Sample ID and VCF file

	output:
	path "*.txt", emit: stats  // Statistics in text format

	script:
	outname = "${vcf.baseName}".replace(".vcf.gz", "")
	"""
	# Skip mitochondrial variants (MT), generate stats for others
	if [[ ! "${vcf}" =~ "_MT" ]]; then
		bcftools stats --threads ${task.cpus} \
		${vcf} > ${outname}.txt
	
	else
		echo "Skipping ${sampleID} MT file." > ${outname}.txt
	fi
	"""
}