/**
 * Minimap2 processes for aligning long reads to a reference genome
 * Includes both indexing the reference and performing the alignment
 */

/**
 * MINIMAP2_INDEX process to create an index of the reference genome
 * This index accelerates subsequent alignments to the reference
 */
process MINIMAP2_INDEX {
	tag "$reference.simpleName"
	label 'process_high'
	debug true

	//publishDir "$params.outdir/minimap2/reference", mode: 'copy'

	input:
	path reference  // Reference genome FASTA

	output:
	path "${filename}.mmi", emit: indexed_reference  // Minimap2 index file
	

	script:
	filename = reference.getSimpleName()  // Get filename without extension
	unzip_reference = reference.toString().replace(".gz", "")  // Remove .gz extension if present
	
	if (reference.toString().endsWith(".gz")) {
		"""
		# For compressed references, uncompress first then index
		gunzip -c $reference > $unzip_reference
		minimap2 -d ${filename}.mmi $unzip_reference
		"""
	} else {
		"""
		# Index uncompressed reference directly
		minimap2 -d ${filename}.mmi $reference
		"""
	}
}

/**
 * MINIMAP2 process to align long reads to a reference genome
 * Optimized for Oxford Nanopore reads with the map-ont preset
 */
process MINIMAP2 {
	tag "$sampleID"
	label "process_high"

	publishDir "$params.outdir/minimap2", mode: 'copy'

	input:
	tuple val( sampleID ), path( reads )  // Sample ID and reads (FASTQ or BAM)
	path reference                        // Reference genome or index

	output:
	tuple val( sampleID ), path( "${sampleID}.sorted.bam" ), path( "${sampleID}.sorted.bam.bai" ), emit: mapped_tuple  // Aligned reads (BAM) with index

	script:
	if (reads.toString().endsWith(".bam")) {
		"""
		# For BAM input, use dorado aligner
		dorado aligner -t $task.cpus --mm2-preset "map-ont" $reference $reads | samtools sort -@ $task.cpus -o ${sampleID}.sorted.bam

		# Index the sorted BAM
		samtools index ${sampleID}.sorted.bam
		"""
	} else {
		"""
		# For FASTQ input, use minimap2 directly
		minimap2 -ax map-ont --MD -t $task.cpus $reference $reads | samtools sort -@ $task.cpus -o ${sampleID}.sorted.bam

		# Index the sorted BAM
		samtools index ${sampleID}.sorted.bam
		"""
	}
}