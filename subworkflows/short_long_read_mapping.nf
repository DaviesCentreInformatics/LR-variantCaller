nextflow.enable.dsl = 2

include { MINIMAP2 as MINIMAP2_MAP            } from '../modules/minimap2'
include { SAMTOOLS_FLAGSTAT 				  } from '../modules/samtools'
include { SAMTOOLS_IDXSTATS                   } from '../modules/samtools'
include { SAMTOOLS_STATS                      } from '../modules/samtools'
include { SAMTOOLS_SPLITBAM as SPLITBAM       } from '../modules/samtools'

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
