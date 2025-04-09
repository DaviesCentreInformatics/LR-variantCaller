/*
 * Workflow: LONG_READ_PREPROCESSING
 * 
 * Purpose:
 *   This workflow performs quality control and filtering of long read data before alignment.
 *   It generates QC reports for both raw and filtered reads using NanoPlot.
 *
 * Inputs:
 *   - samples: Channel containing tuples of sample ID and fastq files
 *
 * Outputs:
 *   - raw_report: QC reports for raw reads from NanoPlot
 *   - filtered_reads: Reads after quality filtering using Filtlong
 *   - filtered_report: QC reports for the filtered reads from NanoPlot
 *
 * Process:
 *   1. Run NanoPlot on raw input reads for QC
 *   2. Filter reads using Filtlong to remove low quality reads
 *   3. Run NanoPlot again on filtered reads to assess filtering effect
 */

nextflow.enable.dsl = 2

include { FILTLONG                      } from '../modules/filtlong' 
include { NANOPLOT as NANOPLOT_RAW      } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_FILTERED } from '../modules/nanoplot'

workflow LONG_READ_PREPROCESSING {
	take:
		
		samples

	main:
		// Generate QC report on raw reads
		NANOPLOT_RAW(samples)
		// Filter reads to remove low quality sequences
		FILTLONG(samples)
		// Generate QC report on filtered reads
		NANOPLOT_FILTERED(FILTLONG.out.result_tuple)
		filtered_report = NANOPLOT_FILTERED.out.report
		filtered_reads  = FILTLONG.out.result_tuple
		
	emit:
		// Output both raw and filtered reports to track effect of filtering
		raw_report      = NANOPLOT_RAW.out.report
		filtered_reads
		filtered_report
}