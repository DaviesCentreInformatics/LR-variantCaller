nextflow.enable.dsl=2

include { SNIFFLES2                     } from '../modules/sniffles2'
include { SVIM                          } from '../modules/svim'
include { CUTESV                        } from '../modules/cutesv'
include { DYSGU                         } from '../modules/dysgu'



workflow LONG_READ_SV_CALLING {
	take:
		mapped
		reference_genome
		reference_genome_index

	main:

		//bam = SPLITBAM(mapped).split_bam.transpose()
		bam = mapped
		// Call SVs
		if (params.sniffles) {
			sniffles = SNIFFLES2(bam, reference_genome, reference_genome_index).res_tuple
		}
		else {
			sniffles = Channel.empty()
		}
		if (params.svim) {
			svim = SVIM(bam, reference_genome, reference_genome_index).res_tuple
		}
		else {
			svim = Channel.empty()
		}
		if (params.cutesv) {
			cutesv = CUTESV(bam, reference_genome, reference_genome_index).res_tuple
		}
		else {
			cutesv = Channel.empty()
		}
		if (params.dysgu) {
			dysgu = DYSGU(bam, reference_genome, reference_genome_index).res_tuple
		}
		else {
			dysgu = Channel.empty()
		}
		
	
	emit:
		sniffles
		svim
		cutesv
		dysgu
}