nextflow.enable.dsl = 2

include { MULTIQC                       } from '../modules/multiqc'
include { SAMTOOLS_FAIDX                } from '../modules/samtools'
include { LONG_READ_VARIANT_CALLING     } from '../subworkflows/long_read_variant_calling'
include { BCFTOOLS_STATS as SNP_STATS   } from '../modules/bcftools'
include { BCFTOOLS_STATS as SNIF_STATS  } from '../modules/bcftools'
include { BCFTOOLS_STATS as SVIM_STATS  } from '../modules/bcftools'
include { BCFTOOLS_STATS as CUSV_STATS  } from '../modules/bcftools'
include { BCFTOOLS_STATS as DYSGU_STATS } from '../modules/bcftools'

/*
 * DEFINE THE MAIN WORKFLOW
 */
workflow ALREADY_MAPPED {
	take:
		samples
		
	main:	
		
		bam = samples

		fasta = params.reference
		fai = params.reference_idx

		LONG_READ_VARIANT_CALLING(bam, fasta, fai)
	
		// snp_stats    = SNP_STATS(LONG_READ_VARIANT_CALLING.out.snps)
		snif_stats   = SNIF_STATS(LONG_READ_VARIANT_CALLING.out.sniffles)
		svim_stats   = SVIM_STATS(LONG_READ_VARIANT_CALLING.out.svim)
		cutesv_stats = CUSV_STATS(LONG_READ_VARIANT_CALLING.out.cutesv)
		dysgu_stats  = DYSGU_STATS(LONG_READ_VARIANT_CALLING.out.dysgu)

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
    
    emit:
        snif_stats
        svim_stats
        cutesv_stats
        dysgu_stats
}