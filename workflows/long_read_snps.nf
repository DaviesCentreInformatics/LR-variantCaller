nextflow.enable.dsl = 2

include { LONG_READ_PREPROCESSING   } from '../subworkflows/long_read_preprocessing'
include { LONG_READ_MAPPING         } from '../subworkflows/long_read_mapping'
include { MULTIQC                   } from '../modules/multiqc'
include { SAMTOOLS_FAIDX            } from '../modules/samtools'
include { LONG_READ_VARIANT_CALLING } from '../subworkflows/modified_long_read_variant_calling'

/*
 * DEFINE THE MAIN WORKFLOW
 */
workflow LONG_READ_VARIANTS {
	take:
		samples
		fasta
		fai
		
	main:	
		
		LONG_READ_VARIANT_CALLING(samples, fasta, fai)
}