/*
 * Workflow for variant calling from pre-mapped BAM files
 * This workflow skips preprocessing and mapping steps,
 * assuming BAM files are already available and properly processed.
 */
nextflow.enable.dsl = 2

// Import required modules and subworkflows
include { MULTIQC                       } from '../modules/multiqc'                        // QC reporting
include { SAMTOOLS_FAIDX                } from '../modules/samtools'                       // Reference indexing
include { LONG_READ_VARIANT_CALLING     } from '../subworkflows/long_read_variant_calling' // Variant calling pipeline
// Import bcftools for different variant callers' stats
include { BCFTOOLS_STATS as SNP_STATS   } from '../modules/bcftools'  // For small variant statistics
include { BCFTOOLS_STATS as SNIF_STATS  } from '../modules/bcftools'  // For Sniffles SV statistics
include { BCFTOOLS_STATS as SVIM_STATS  } from '../modules/bcftools'  // For SVIM SV statistics
include { BCFTOOLS_STATS as CUSV_STATS  } from '../modules/bcftools'  // For CuteSV SV statistics
include { BCFTOOLS_STATS as DYSGU_STATS } from '../modules/bcftools'  // For Dysgu SV statistics

/*
 * DEFINE THE MAIN WORKFLOW
 */
workflow ALREADY_MAPPED {
	take:
		samples  // Input channel with BAM files
		
	main:	
		// Use the provided BAM files directly
		bam = samples

		// Use provided reference files
		fasta = params.reference       // Reference genome fasta
		fai = params.reference_idx     // Reference genome index

		// Perform variant calling using multiple variant callers
		LONG_READ_VARIANT_CALLING(bam, fasta, fai)
	
		// Generate statistics for the variant calls
		// SNP statistics currently disabled
		// snp_stats    = SNP_STATS(LONG_READ_VARIANT_CALLING.out.snps)
		
		// Generate statistics for structural variants
		snif_stats   = SNIF_STATS(LONG_READ_VARIANT_CALLING.out.sniffles)  // Sniffles SV statistics
		svim_stats   = SVIM_STATS(LONG_READ_VARIANT_CALLING.out.svim)      // SVIM SV statistics
		cutesv_stats = CUSV_STATS(LONG_READ_VARIANT_CALLING.out.cutesv)    // CuteSV SV statistics
		dysgu_stats  = DYSGU_STATS(LONG_READ_VARIANT_CALLING.out.dysgu)    // Dysgu SV statistics

		// MultiQC report generation is disabled since we don't have preprocessing/mapping metrics
		// reports_and_logs = Channel.empty()
		// //reports_and_logs.view()
		// multiqc_input_ch = reports_and_logs.mix(
		// 	LONG_READ_PREPROCESSING.out.raw_report,
		// 	LONG_READ_PREPROCESSING.out.filtered_report,
		// 	LONG_READ_MAPPING.out.stats,
		// 	LONG_READ_MAPPING.out.idxstat,
		// 	LONG_READ_MAPPING.out.flagstat,
		// 	LONG_READ_MAPPING.out.coverage)
		// 	.collect()
		// //multiqc_input_ch.view()
		// MULTIQC(multiqc_input_ch)
    
    // Emit outputs for potential downstream processes
    emit:
        snif_stats   // Sniffles statistics
        svim_stats   // SVIM statistics 
        cutesv_stats // CuteSV statistics
        dysgu_stats  // Dysgu statistics
}