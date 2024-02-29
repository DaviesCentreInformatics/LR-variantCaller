nextflow.enable.dsl=2

include { SNIFFLES2               } from '../modules/sniffles2'


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