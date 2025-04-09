/*
 * Workflow: LONG_READ_MAPPING (Archived)
 * 
 * Purpose:
 *   This is an archived version of the mapping workflow that maps long reads to a reference genome
 *   using Minimap2 and processes the alignments with Samtools.
 *
 * Inputs:
 *   - filtered_reads: Channel containing filtered long read data
 *   - minimap_index: Pre-built Minimap2 index
 *
 * Outputs:
 *   - mapped_reads: BAM files split by chromosome for parallel processing
 *   - flagstat: Mapping statistics from Samtools flagstat
 *   - idxstat: Index statistics from Samtools idxstats
 *   - stats: Detailed alignment statistics from Samtools stats
 * 
 * Note:
 *   This version always splits BAM files, which is now optional in newer workflows.
 */

nextflow.enable.dsl = 2

include { MINIMAP2 as MINIMAP2_MAP            } from '../../modules/minimap2'
include { SAMTOOLS_FLAGSTAT 				  } from '../../modules/samtools'
include { SAMTOOLS_IDXSTATS                   } from '../../modules/samtools'
include { SAMTOOLS_STATS                      } from '../../modules/samtools'
include { SAMTOOLS_SPLITBAM as SPLITBAM       } from '../../modules/samtools'

workflow LONG_READ_MAPPING {
	take:
	filtered_reads
	minimap_index

	main:
		MINIMAP2_MAP(filtered_reads, minimap_index)
		mapped = MINIMAP2_MAP.out.mapped_tuple
		
		SAMTOOLS_FLAGSTAT(mapped)
		SAMTOOLS_IDXSTATS(mapped)
		SAMTOOLS_STATS(mapped)
		SPLITBAM(mapped)
		
	emit:
		mapped_reads = SPLITBAM.out.split_bam.transpose()
		flagstat = SAMTOOLS_FLAGSTAT.out
		idxstat = SAMTOOLS_IDXSTATS.out
		stats = SAMTOOLS_STATS.out
}
