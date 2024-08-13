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
		
	main:	
		LONG_READ_PREPROCESSING(samples)

		bam = LONG_READ_MAPPING(LONG_READ_PREPROCESSING.out.filtered_reads,
	                  params.minimap_index).mapped
		
		fasta = params.reference
		fai = params.reference_idx

		LONG_READ_VARIANT_CALLING(samples, fasta, fai)
}