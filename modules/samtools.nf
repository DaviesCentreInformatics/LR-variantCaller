/**
 * Collection of Samtools processes for BAM file manipulation and analysis
 * Includes procedures for BAM statistics, indexing, flagstat, and more
 */

/**
 * SAMTOOLS_FLAGSTAT process to generate alignment statistics 
 * Counts the number of alignments for each FLAG type
 */
process SAMTOOLS_FLAGSTAT {
	tag "$sampleID"
	label "process_low"

	publishDir "$params.outdir/samtools", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("flagstats") > 0 ) "flagstats/$filename"
			else null }

	input:
	tuple val(sampleID), path(bam), path(bai)  // BAM file with index

	output:
	path "${sampleID}.flagstats", emit: flagstats  // Flagstat output file


	script:
	"""
	# Generate alignment statistics using flagstat
	samtools flagstat -@ ${task.cpus} ${bam} > ${sampleID}.flagstats
	"""
}

/**
 * SAMTOOLS_IDXSTATS process to retrieve statistics from the BAM index
 * Reports counts of alignments per reference sequence
 */
process SAMTOOLS_IDXSTATS {
	tag "$sampleID"
	label "process_low"

	publishDir "$params.outdir/samtools", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("idxstats") > 0 ) "idxstats/$filename"
			else null }

	input:
	tuple val(sampleID), path(bam), path(bai)  // BAM file with index

	output:
	path "${sampleID}.idxstats", emit: idxstats  // Idxstats output file


	script:
	"""
	# Generate statistics from BAM index
	samtools idxstats -@ ${task.cpus} ${bam} > ${sampleID}.idxstats
	"""
}

/**
 * SAMTOOLS_STATS process to generate comprehensive statistics for a BAM file
 * Provides detailed information about alignments, mapping quality, etc.
 */
process SAMTOOLS_STATS {
	tag "$sampleID"
	label "process_low"

	publishDir "$params.outdir/samtools", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("stats") > 0 ) "stats/$filename"
			else null }

	input:
	tuple val(sampleID), path(bam), path(bai)  // BAM file with index

	output:
	path "${sampleID}.stats", emit: stats  // Stats output file


	script:
	"""
	# Generate comprehensive BAM statistics
	samtools stats -@ ${task.cpus} ${bam} > ${sampleID}.stats
	"""
}

/**
 * SAMTOOLS_SPLITBAM process to split a BAM file by chromosome
 * Creates separate BAM files for each main chromosome
 */
process SAMTOOLS_SPLITBAM {
	tag "$sampleID"
	label "process_low"

	//publishDir "$params.outdir/variants/splitBams", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)  // BAM file with index

	output:
    tuple val(sampleID), path("${sampleID}_*.bam"), path("${sampleID}_*.bam.bai"), emit: split_bam  // Split BAMs with indexes

	shell:
    '''
    # Extract standard chromosomes from the BAM index (1-22, X, Y)
    samtools idxstats -@ !{task.cpus} !{bam} | cut -f 1 | grep -E '^([0-9]{1,2}|X|Y)$' > !{sampleID}.chromosomes.txt
    
    # For each chromosome, extract reads and create an indexed BAM file
    while IFS= read -r line; do
        samtools view -@ !{task.cpus} -b !{bam} ${line} > !{sampleID}_${line}.bam ;
        samtools index -@ !{task.cpus} !{sampleID}_${line}.bam
    done < !{sampleID}.chromosomes.txt
    '''
}

/**
 * SAMTOOLS_FAIDX process to index a FASTA reference genome
 * Creates an index (.fai) for efficient random access to the reference
 */
process SAMTOOLS_FAIDX {
	tag "$reference.simpleName"
	label "process_low"

	input:
	path reference  // FASTA reference genome (possibly compressed)

	output:
	path "${unzip_reference}"          , emit: fa   // Uncompressed FASTA
	path "${unzip_reference}.fai"      , emit: fai  // FASTA index

	script:
	unzip_reference = reference.toString().replace(".gz", "")  // Remove .gz extension if present

	if (reference.toString().endsWith(".gz")) {
		"""
		# For compressed references, uncompress first then index
		gunzip -c ${reference} > ${unzip_reference}
		samtools faidx ${unzip_reference}
		"""
	} else {
		"""
		# Index uncompressed reference directly
		samtools faidx ${reference}
		"""
	}
}

/**
 * SAMTOOLS_INDEX process to create an index for a BAM file
 * Enables random access to the BAM file for efficient analysis
 */
process SAMTOOLS_INDEX {
	tag "$sampleID"
	label "process_low"

	input:
	tuple val(sampleID), path(bam)  // BAM file without index

	output:
	tuple val(sampleID), path(bam), path(bam + ".bai"), emit: bam  // BAM with its index


	script:
	"""
	# Create BAM index
	samtools index -@ ${task.cpus} ${bam}
	"""
}