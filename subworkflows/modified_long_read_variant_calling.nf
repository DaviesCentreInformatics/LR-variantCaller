nextflow.enable.dsl=2

include { CLAIR3                        } from '../modules/clair3'
include { SAMTOOLS_SPLITBAM as SPLITBAM } from '../modules/samtools'



workflow LONG_READ_VARIANT_CALLING {
	take:
		bam
		reference_genome
		reference_genome_index

	main:
		// Split the BAM
		split_bams = SPLITBAM(bam).split_bam.transpose()
		// Call SNPs
		CLAIR3(split_bams, reference_genome, reference_genome_index)

		// Call SVs
		//SNIFFLES2(bam, reference_genome, reference_genome_index)
	
	emit:
		// snps = CLAIR3.out
		// svs  = SNIFFLES2.out
}