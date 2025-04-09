/*
 * Workflow for SNP/small variant calling from long reads
 * This workflow handles preprocessing, mapping and small variant calling
 * but excludes structural variant calling.
 */
nextflow.enable.dsl = 2

// Import required subworkflows and modules
include { LONG_READ_PREPROCESSING   } from '../subworkflows/long_read_preprocessing'    // QC and filtering
include { LONG_READ_MAPPING         } from '../subworkflows/long_read_mapping'          // Read alignment
include { MULTIQC                   } from '../modules/multiqc'                         // QC reporting
include { SAMTOOLS_FAIDX            } from '../modules/samtools'                        // Reference indexing
include { LONG_READ_VARIANT_CALLING } from '../subworkflows/archive/modified_long_read_variant_calling'  // Using modified variant calling

/*
 * DEFINE THE MAIN WORKFLOW
 */
workflow LONG_READ_VARIANTS {
	take:
		samples  // Input channel of samples containing fastq files
		
	main:	
		// Step 1: Preprocess the long reads (QC and filtering)
		LONG_READ_PREPROCESSING(samples)

		// Step 2: Map filtered reads to reference genome
		bam = LONG_READ_MAPPING(LONG_READ_PREPROCESSING.out.filtered_reads,
	                  params.minimap_index).mapped
		
		// Use provided reference files
		fasta = params.reference       // Reference genome fasta
		fai = params.reference_idx     // Reference genome index

		// Step 3: Call small variants using the modified variant calling pipeline
		// This uses a customized pipeline for SNPs and small indels only
		LONG_READ_VARIANT_CALLING(bam, fasta, fai)
}