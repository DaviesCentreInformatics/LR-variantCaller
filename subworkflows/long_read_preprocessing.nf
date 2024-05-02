nextflow.enable.dsl = 2

include { FILTLONG                      } from '../modules/filtlong' 
include { NANOPLOT as NANOPLOT_RAW      } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_FILTERED } from '../modules/nanoplot'
include { FILTER_READS                  } from '../modules/filter_reads'

workflow LONG_READ_PREPROCESSING {
	take: samples

	main:
		NANOPLOT_RAW(samples)
		if (!params.skip_filtlong) {
			FILTLONG(samples)
			NANOPLOT_FILTERED(FILTLONG.out.result_tuple)
			filtered_report = NANOPLOT_FILTERED.out.report
			//filtered_reads  = FILTLONG.out.result_tuple
			filtered_reads  = FILTER_READS(filtered_reads)
		}
		else {
			filtered_reads = samples
			filtered_report = Channel.empty()
		}
	
	emit:
		raw_report      = NANOPLOT_RAW.out.report
		filtered_reads
		filtered_report
}