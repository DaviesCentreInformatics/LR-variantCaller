/*
 * Main workflow for long read variant calling
 * This workflow handles the complete process from preprocessing raw reads
 * through mapping and variant calling.
 */
nextflow.enable.dsl = 2

// Import required modules and subworkflows
include { LONG_READ_PREPROCESSING       } from '../subworkflows/long_read_preprocessing'  // QC and filtering of raw reads
include { LONG_READ_MAPPING             } from '../subworkflows/long_read_mapping'        // Alignment to reference genome
include { MULTIQC                       } from '../modules/multiqc'                       // Quality control report generation
include { SAMTOOLS_FAIDX                } from '../modules/samtools'                      // Reference genome indexing
include { LONG_READ_VARIANT_CALLING     } from '../subworkflows/long_read_variant_calling'// Variant calling using multiple callers
// Import bcftools for different variant callers' stats
include { BCFTOOLS_STATS as SNP_STATS   } from '../modules/bcftools'  // For small variant statistics
include { BCFTOOLS_STATS as SNIF_STATS  } from '../modules/bcftools'  // For Sniffles SV statistics
include { BCFTOOLS_STATS as SVIM_STATS  } from '../modules/bcftools'  // For SVIM SV statistics
include { BCFTOOLS_STATS as CUSV_STATS  } from '../modules/bcftools'  // For CuteSV SV statistics
include { BCFTOOLS_STATS as DYSGU_STATS } from '../modules/bcftools'  // For Dysgu SV statistics


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
		//bams.view()  // Debugging line (commented out)
		
		// Use provided reference files instead of generating indices
		// (fasta, fai) = SAMTOOLS_FAIDX(params.reference)  // Commented out as using pre-indexed reference
		fasta = params.reference       // Reference genome fasta
		fai = params.reference_idx     // Reference genome index

		// Step 3: Call variants using multiple variant callers
		LONG_READ_VARIANT_CALLING(bam, fasta, fai)
	
		// Generate statistics for the variant calls
		// snp_stats    = SNP_STATS(LONG_READ_VARIANT_CALLING.out.snps)  // Currently disabled
		snif_stats   = SNIF_STATS(LONG_READ_VARIANT_CALLING.out.sniffles)  // Sniffles SV statistics
		svim_stats   = SVIM_STATS(LONG_READ_VARIANT_CALLING.out.svim)      // SVIM SV statistics
		cutesv_stats = CUSV_STATS(LONG_READ_VARIANT_CALLING.out.cutesv)    // CuteSV SV statistics
		dysgu_stats  = DYSGU_STATS(LONG_READ_VARIANT_CALLING.out.dysgu)    // Dysgu SV statistics

		// Step 4: Prepare and generate MultiQC report from all QC outputs
		reports_and_logs = Channel.empty()  // Initialize empty channel for QC reports
		//reports_and_logs.view()  // Debugging line (commented out)
		
		// Combine all QC reports for MultiQC
		multiqc_input_ch = reports_and_logs.mix(
			LONG_READ_PREPROCESSING.out.raw_report,      // Raw read QC metrics
			LONG_READ_PREPROCESSING.out.filtered_report, // Filtered read QC metrics
			LONG_READ_MAPPING.out.stats,                 // Mapping statistics
			LONG_READ_MAPPING.out.idxstat,               // Index statistics
			LONG_READ_MAPPING.out.flagstat,              // Flag statistics
			LONG_READ_MAPPING.out.coverage)              // Coverage metrics
			.collect()
		//multiqc_input_ch.view()  // Debugging line (commented out)
		
		// Generate comprehensive QC report
		MULTIQC(multiqc_input_ch)
}