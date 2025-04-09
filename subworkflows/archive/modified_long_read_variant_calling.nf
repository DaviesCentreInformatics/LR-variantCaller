/*
 * Workflow: LONG_READ_VARIANT_CALLING (Archived)
 * 
 * Purpose:
 *   This is an archived version of the variant calling workflow that focuses on
 *   SNP/small variant calling using Clair3 after splitting BAM files by chromosome.
 *
 * Inputs:
 *   - bam: Channel containing BAM files with long reads aligned to reference
 *   - reference_genome: Reference genome FASTA file
 *   - reference_genome_index: Index file for the reference genome
 *
 * Process:
 *   1. Split BAM files by chromosome for parallel processing
 *   2. Run Clair3 on each chromosome for SNP/small variant calling
 *
 * Note:
 *   This version has SV calling commented out and doesn't emit outputs properly.
 *   It has been superseded by more comprehensive workflows.
 */

nextflow.enable.dsl=2

include { CLAIR3                        } from '../modules/clair3'
include { SAMTOOLS_SPLITBAM as SPLITBAM } from '../modules/samtools'



workflow LONG_READ_VARIANT_CALLING {
	take:
		bam
		reference_genome
		reference_genome_index

	main:
		// Split the BAM file by chromosome for parallel processing
		split_bams = SPLITBAM(bam).split_bam.transpose()
		// Call SNPs with Clair3
		CLAIR3(split_bams, reference_genome, reference_genome_index)

		// Call SVs
		//SNIFFLES2(bam, reference_genome, reference_genome_index)
	
	// emit:
	// 	// snps = CLAIR3.out
	// 	// svs  = SNIFFLES2.out
}