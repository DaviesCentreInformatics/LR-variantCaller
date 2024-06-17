process SURVIVOR {
	tag "$sampleID"
	label "process_low","error_retry"

	debug true

	publishDir "$params.outdir/variants/SVs/SURVIVOR", mode: 'copy'

	input:
	tuple val(sampleID), path(vcfs)

	output:
	path "${sampleID}_survivor.vcf.gz", emit: merged_svs
	path "*.txt",                       emit: survivor_input
	//stdout

	script:
	vcf_files = vcfs.collect { "${it}"}
	"""
	createSurvivorInput.py --sampleID ${sampleID} \
		--output surivor_input.txt \
		${vcf_files.join(" ")}
	survivor merge ${sampleID}_survivor_input.txt 10 2 1 0 0 50 \
		${sampleID}_survivor.vcf
	pigz ${sampleID}_survivor.vcf
	"""
}