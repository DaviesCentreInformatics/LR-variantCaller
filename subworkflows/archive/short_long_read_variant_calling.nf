/*
 * Workflow: LONG_READ_VARIANT_CALLING (Archived)
 * 
 * Purpose:
 *   This is an archived version of the variant calling workflow that only uses Sniffles2
 *   for structural variant detection from long reads.
 *
 * Inputs:
 *   - bam: Channel containing BAM files
 *   - reference_genome: Reference genome FASTA file
 *   - reference_genome_index: Index file for the reference genome
 *
 * Outputs:
 *   - svs: Structural variants called by Sniffles2
 *
 * Note:
 *   This is a simplified early version that has been superseded by more comprehensive workflows.
 */

nextflow.enable.dsl=2

include { SNIFFLES2               } from '../../modules/sniffles2'


workflow LONG_READ_VARIANT_CALLING {
	take:
		bam
		reference_genome
		reference_genome_index

	main:
		SNIFFLES2(bam, reference_genome, reference_genome_index)
	
	emit:
		svs = SNIFFLES2.out
}