nextflow.enable.dsl=2

include { CLAIR3                      } from '../modules/clair3'
include { SNIFFLES2                   } from '../modules/sniffles2'
include { BCFTOOLS_STATS as SV_STATS  } from '../modules/bcftools'
include { BCFTOOLS_STATS as SNP_STATS } from '../modules/bcftools'


workflow LONG_READ_VARIANT_CALLING {
	take:
		bam
		reference_genome
		reference_genome_index

	main:
		CLAIR3(bam, reference_genome, reference_genome_index)
		SNIFFLES2(bam, reference_genome, reference_genome_index)
		
		SV_STATS(SNIFFLES2.out.sv_vcf, reference_genome, reference_genome_index)
	
	emit:
		snps = CLAIR3.out
		svs  = SNIFFLES2.out[1]
		sv_stats = SV_STATS.out
}