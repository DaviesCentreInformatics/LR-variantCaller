nextflow.enable.dsl = 2

include { LONG_READ_PREPROCESSING   } from '../subworkflows/long_read_preprocessing'
include { LONG_READ_MAPPING         } from '../subworkflows/long_read_mapping'
include { MULTIQC                   } from '../modules/multiqc'
include { SAMTOOLS_FAIDX            } from '../modules/samtools'
include { LONG_READ_VARIANT_CALLING } from '../subworkflows/long_read_variant_calling'

/*
 * DEFINE THE MAIN WORKFLOW
 */
workflow LONG_READ_VARIANTS {
	take:
		samples
		
	main:	
		LONG_READ_PREPROCESSING(samples)
	
		bams = LONG_READ_MAPPING(LONG_READ_PREPROCESSING.out.filtered_reads,
	                  params.minimap_index).mapped_reads
		//bams.view()
		(fasta, fai) = SAMTOOLS_FAIDX(params.reference)

		LONG_READ_VARIANT_CALLING(bams, fasta, fai)
	
		reports_and_logs = Channel.empty()
		//reports_and_logs.view()
		multiqc_input_ch = reports_and_logs.mix(LONG_READ_PREPROCESSING.out.raw_report,
	                                        	LONG_READ_PREPROCESSING.out.filtered_report,
												LONG_READ_MAPPING.out.stats,
												LONG_READ_MAPPING.out.idxstat,
												LONG_READ_MAPPING.out.flagstat,
												LONG_READ_MAPPING.out.coverage,
												LONG_READ_VARIANT_CALLING.out.snps,
												LONG_READ_VARIANT_CALLING.out.svs)
						                	.collect()
		//multiqc_input_ch.view()
		MULTIQC(multiqc_input_ch)
	
	
}