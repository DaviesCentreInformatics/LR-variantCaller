process SAMTOOLS_FLAGSTAT {
	tag "$sampleID"
	label "process_low"

	publishDir "$params.outdir/samtools", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("flagstats") > 0 ) "flagstats/$filename"
			else null }

	input:
	tuple val(sampleID), path(bam), path(bai)

	output:
	path "${sampleID}.flagstats", emit: flagstats


	script:
	"""
	samtools flagstat -@ ${task.cpus} ${bam} > ${sampleID}.flagstats
	"""
}

process SAMTOOLS_IDXSTATS {
	tag "$sampleID"
	label "process_low"

	publishDir "$params.outdir/samtools", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("idxstats") > 0 ) "idxstats/$filename"
			else null }

	input:
	tuple val(sampleID), path(bam), path(bai)

	output:
	path "${sampleID}.idxstats", emit: idxstats


	script:
	"""
	samtools idxstats -@ ${task.cpus} ${bam} > ${sampleID}.idxstats
	"""
}

process SAMTOOLS_STATS {
	tag "$sampleID"
	label "process_low"

	publishDir "$params.outdir/samtools", mode: 'copy',
		saveAs: { filename -> 
			if (filename.indexOf("stats") > 0 ) "stats/$filename"
			else null }

	input:
	tuple val(sampleID), path(bam), path(bai)

	output:
	path "${sampleID}.stats", emit: stats


	script:
	"""
	samtools stats -@ ${task.cpus} ${bam} > ${sampleID}.stats
	"""
}

process SAMTOOLS_SPLITBAM {
	tag "$sampleID"
	label "process_low"

	//publishDir "$params.outdir/variants/splitBams", mode: 'copy'

	input:
	tuple val(sampleID), path(bam), path(bai)

	output:
    tuple val(sampleID), path("${sampleID}_*.bam"), path("${sampleID}_*.bam.bai"), emit: split_bam

	shell:
    '''
    samtools idxstats -@ !{task.cpus} !{bam} | cut -f 1 | grep -E '^([0-9]{1,2}|X|Y)$' > !{sampleID}.chromosomes.txt
    while IFS= read -r line; do
        samtools view -@ !{task.cpus} -b !{bam} ${line} > !{sampleID}_${line}.bam ;
        samtools index -@ !{task.cpus} !{sampleID}_${line}.bam
    done < !{sampleID}.chromosomes.txt
    '''
}

process SAMTOOLS_FAIDX {
	tag "$reference.simpleName"
	label "process_low"

	input:
	path reference

	output:
	path "${unzip_reference}"          , emit: fa
	path "${unzip_reference}.fai"      , emit: fai

	script:
	unzip_reference = reference.toString().replace(".gz", "")

	if (reference.toString().endsWith(".gz")) {
		"""
		gunzip -c ${reference} > ${unzip_reference}
		samtools faidx ${unzip_reference}
		"""
	} else {
		"""
		samtools faidx ${reference}
		"""
	}
}
