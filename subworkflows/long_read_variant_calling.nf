nextflow.enable.dsl=2

include { CLAIR3                  } from '../modules/clair3'
include { SNIFFLES2               } from '../modules/sniffles2'
include { SVIM                    } from '../modules/svim'
include { CUTESV                  } from '../modules/cutesv'
include { DYSGU                   } from '../modules/dysgu'



workflow LONG_READ_VARIANT_CALLING {
	take:
		bam
		reference_genome
		reference_genome_index

	main:
		// Call SNPs
		if (params.skip_snps == false) {
			CLAIR3(bam, reference_genome, reference_genome_index)
			snps = CLAIR3.out.gvcf
		} else {
			snps = Channel.empty()
		}

		// Call SVs
		SNIFFLES2(bam, reference_genome, reference_genome_index)
		SVIM(bam, reference_genome, reference_genome_index)
		CUTESV(bam, reference_genome, reference_genome_index)
		DYSGU(bam, reference_genome, reference_genome_index)
	
	emit:
		snps
		sniffles  = SNIFFLES2.out.res_tuple
		svim = SVIM.out.res_tuple
		cutesv = CUTESV.out.res_tuple
		dysgu = DYSGU.out.res_tuple
}