nextflow.enable.dsl = 2

include { LONG_READ_PREPROCESSING       } from '../subworkflows/long_read_preprocessing'
include { LONG_READ_MAPPING             } from '../subworkflows/long_read_mapping'
include { MULTIQC                       } from '../modules/multiqc'
include { SAMTOOLS_FAIDX                } from '../modules/samtools'
include { LONG_READ_VARIANT_CALLING     } from '../subworkflows/long_read_variant_calling'
include { BCFTOOLS_STATS as SNP_STATS   } from '../modules/bcftools'
include { BCFTOOLS_STATS as SNIF_STATS  } from '../modules/bcftools'
include { BCFTOOLS_STATS as SVIM_STATS  } from '../modules/bcftools'
include { BCFTOOLS_STATS as CUSV_STATS  } from '../modules/bcftools'
include { BCFTOOLS_STATS as DYSGU_STATS } from '../modules/bcftools'
include { MODKIT						} from '../modules/modkit'

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
		//bams.view()
		// (fasta, fai) = SAMTOOLS_FAIDX(params.reference)
		fasta = params.reference
		fai = params.reference_idx

		LONG_READ_VARIANT_CALLING(bam, fasta, fai)

		//svs_to_merge = Channel.empty()
		// svs_to_merge = Channel.empty()
		
		// survivor_input = svs_to_merge.mix(LONG_READ_VARIANT_CALLING.out.sniffles,
		// 										LONG_READ_VARIANT_CALLING.out.svim,
		// 										LONG_READ_VARIANT_CALLING.out.cutesv,
		// 										LONG_READ_VARIANT_CALLING.out.dysgu).groupTuple(by: 0)
		
	
		// snp_stats    = SNP_STATS(LONG_READ_VARIANT_CALLING.out.snps)
		snif_stats   = SNIF_STATS(LONG_READ_VARIANT_CALLING.out.sniffles)
		svim_stats   = SVIM_STATS(LONG_READ_VARIANT_CALLING.out.svim)
		cutesv_stats = CUSV_STATS(LONG_READ_VARIANT_CALLING.out.cutesv)
		dysgu_stats  = DYSGU_STATS(LONG_READ_VARIANT_CALLING.out.dysgu)

		reports_and_logs = Channel.empty()
		//reports_and_logs.view()
		multiqc_input_ch = reports_and_logs.mix(
			LONG_READ_PREPROCESSING.out.raw_report,
			LONG_READ_PREPROCESSING.out.filtered_report,
			LONG_READ_MAPPING.out.stats,
			LONG_READ_MAPPING.out.idxstat,
			LONG_READ_MAPPING.out.flagstat,
			LONG_READ_MAPPING.out.coverage)
			.collect()
		//multiqc_input_ch.view()
		MULTIQC(multiqc_input_ch)
}