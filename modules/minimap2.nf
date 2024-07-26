process MINIMAP2_INDEX {
	tag "$reference.simpleName"
	label 'process_high'
	debug true

	//publishDir "$params.outdir/minimap2/reference", mode: 'copy'

	input:
	path reference

	output:
	path "${filename}.mmi"                               ,  emit: indexed_reference
	

	script:
	filename = reference.getSimpleName()
	unzip_reference = reference.toString().replace(".gz", "")
	
	if (reference.toString().endsWith(".gz")) {
		"""
		gunzip -c $reference > $unzip_reference
		minimap2 -d ${filename}.mmi $unzip_reference
		"""
	} else {
		"""
		minimap2 -d ${filename}.mmi $reference
		"""
	}
}

process MINIMAP2 {
	tag "$sampleID"
	label "process_high"

	publishDir "$params.outdir/minimap2", mode: 'copy'

	input:
	tuple val( sampleID ), path( reads )
	path reference

	output:
	tuple val( sampleID ), path( "${sampleID}.sorted.bam" ), path( "${sampleID}.sorted.bam.bai" ), emit: mapped_tuple

	script:
	if (reads.toString().endsWith(".bam")) {
		"""
		dorado aligner -t $task.cpus --mm2-preset "map-ont" $reference $reads | samtools sort -@ $task.cpus -o ${sampleID}.sorted.bam

		samtools index ${sampleID}.sorted.bam
		"""
	} else {
		"""
		minimap2 -ax map-ont --MD -t $task.cpus $reference $reads | samtools sort -@ $task.cpus -o ${sampleID}.sorted.bam

		samtools index ${sampleID}.sorted.bam
		"""
	}
}