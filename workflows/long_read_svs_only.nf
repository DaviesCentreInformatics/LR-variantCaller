nextflow.enable.dsl = 2

include { LONG_READ_SV_CALLING     } from '../subworkflows/long_read_SV_calling'
include { SAMTOOLS_FAIDX           } from '../modules/samtools'


/*
 * DEFINE THE MAIN WORKFLOW
 */
workflow LONG_READ_VARIANTS {
	take:
		samples
		fasta
		fai

	main:	
		
		//SAMTOOLS_FAIDX(ref)

		LONG_READ_SV_CALLING(bams, fasta, fai)

}