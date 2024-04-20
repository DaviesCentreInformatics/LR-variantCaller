nextflow.enable.dsl=2

include { CLAIR3                  } from '../modules/clair3'
include { SNIFFLES2               } from '../modules/sniffles2'


workflow LONG_READ_VARIANT_CALLING {
	take:
		bam
		reference_genome
		reference_genome_index

	main:
		// Call SNPs
		CLAIR3(bam, reference_genome, reference_genome_index)

		// Call SVs
		SNIFFLES2(bam, reference_genome, reference_genome_index)
	
	emit:
		snps = CLAIR3.out.gvcf
		svs  = SNIFFLES2.out
}