nextflow.enable.dsl = 2

include { FILTLONG                      } from '../modules/filtlong' 
include { NANOPLOT as NANOPLOT_RAW      } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_FILTERED } from '../modules/nanoplot'

workflow LONG_READ_PREPROCESSING {
	take: samples

	main:
		NANOPLOT_RAW(samples)
		FILTLONG(samples)
		NANOPLOT_FILTERED(FILTLONG.out.result_tuple)
	
	emit:
		filtered_reads  = FILTLONG.out.result_tuple
		raw_report      = NANOPLOT_RAW.out.report
		filtered_report = NANOPLOT_FILTERED.out.report
		
}