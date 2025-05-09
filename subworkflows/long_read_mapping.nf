/*
 * Workflow: LONG_READ_MAPPING
 * 
 * Purpose:
 *   This workflow maps long reads to a reference genome using Minimap2 or
 *   Dorado and generates various QC metrics for the mapped reads.
 *
 * Inputs:
 *   - filtered_reads: Channel containing filtered long read data
 *   - reference_genome: Reference genome file
 *
 * Outputs:
 *   - mapped: Aligned BAM files
 *   - coverage: Coverage statistics from Mosdepth
 *   - flagstat: Mapping statistics from Samtools flagstat
 *   - idxstat: Index statistics from Samtools idxstats
 *   - stats: Detailed alignment statistics from Samtools stats
 *
 * Parameters:
 *   - params.minimap_index: Pre-built minimap2 index (optional)
 * 
 * Process:
 *   1. Map reads to reference using Minimap2
 *   2. Calculate coverage with Mosdepth
 *   3. Generate mapping statistics with Samtools
 */

nextflow.enable.dsl = 2

include { MINIMAP2_INDEX                      } from '../modules/minimap2'
include { MINIMAP2 as MINIMAP2_MAP            } from '../modules/minimap2'
include { SAMTOOLS_FLAGSTAT 				  } from '../modules/samtools'
include { SAMTOOLS_IDXSTATS                   } from '../modules/samtools'
include { SAMTOOLS_STATS                      } from '../modules/samtools'
include { SAMTOOLS_SPLITBAM as SPLITBAM       } from '../modules/samtools'
include { MOSDEPTH                            } from '../modules/mosdepth'

workflow LONG_READ_MAPPING {
	take:
		filtered_reads
		reference_genome

	main:
		//MINIMAP2_INDEX(reference_genome)
		//index = MINIMAP2_INDEX.out.indexed_reference
		// Use pre-built index if available
		index = params.minimap_index

		// Map filtered reads to reference
		MINIMAP2_MAP(filtered_reads, index)
		mapped = MINIMAP2_MAP.out.mapped_tuple

		// Calculate coverage statistics
		MOSDEPTH(mapped)
		
		// Generate various mapping statistics
		SAMTOOLS_FLAGSTAT(mapped)
		SAMTOOLS_IDXSTATS(mapped)
		SAMTOOLS_STATS(mapped)
		//SPLITBAM(mapped)
		
	emit:
		mapped
		//mapped_reads = SPLITBAM.out.split_bam.transpose()
		coverage = MOSDEPTH.out
		flagstat = SAMTOOLS_FLAGSTAT.out
		idxstat = SAMTOOLS_IDXSTATS.out
		stats = SAMTOOLS_STATS.out
}
