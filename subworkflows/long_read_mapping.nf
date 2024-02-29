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
		MINIMAP2_INDEX(reference_genome)
		index = MINIMAP2_INDEX.out.indexed_reference

		MINIMAP2_MAP(filtered_reads, index)
		mapped = MINIMAP2_MAP.out.mapped_tuple

		MOSDEPTH(mapped)
		
		SAMTOOLS_FLAGSTAT(mapped)
		SAMTOOLS_IDXSTATS(mapped)
		SAMTOOLS_STATS(mapped)
		SPLITBAM(mapped)
		
	emit:
		mapped_reads = SPLITBAM.out.split_bam.transpose()
		coverage = MOSDEPTH.out
		flagstat = SAMTOOLS_FLAGSTAT.out
		idxstat = SAMTOOLS_IDXSTATS.out
		stats = SAMTOOLS_STATS.out
}
