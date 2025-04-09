/*
 * Workflow for structural variant (SV) calling only
 * This workflow takes pre-mapped BAM files and performs SV calling
 * without preprocessing or mapping steps.
 */
nextflow.enable.dsl = 2

// Import required subworkflows and modules
include { LONG_READ_SV_CALLING     } from '../subworkflows/long_read_SV_calling'  // SV calling subworkflow
include { SAMTOOLS_FAIDX           } from '../modules/samtools'                    // For reference indexing


/*
 * DEFINE THE MAIN WORKFLOW
 */
workflow LONG_READ_VARIANTS {
	take:
		samples  // Input channel with BAM files
		fasta    // Reference genome fasta
		fai      // Reference genome index

	main:	
		// Note: SAMTOOLS_FAIDX is commented out as we're using pre-indexed reference
		//SAMTOOLS_FAIDX(ref)

		// Perform structural variant calling using multiple SV callers
		LONG_READ_SV_CALLING(samples, fasta, fai)  // Calls SVs using specified methods

}