nextflow.enable.dsl=2

include { CLAIR3                  } from '../modules/clair3'
include { SNIFFLES2               } from '../modules/sniffles2'


workflow LONG_READ_VARIANT_CALLING {
	take:
		bam
		reference_genome
		reference_genome_index

	main:
		CLAIR3(bam, reference_genome, reference_genome_index)
		SNIFFLES2(bam, reference_genome, reference_genome_index)
	
	emit:
		snps = CLAIR3.out
		svs  = SNIFFLES2.out
}