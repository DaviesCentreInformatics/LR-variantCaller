nextflow.enable.dsl=2

include { CLAIR3                        } from '../modules/clair3'
include { SNIFFLES2                     } from '../modules/sniffles2'
include { SVIM                          } from '../modules/svim'
include { CUTESV                        } from '../modules/cutesv'
include { DYSGU                         } from '../modules/dysgu'
include { SAMTOOLS_SPLITBAM as SPLITBAM } from '../modules/samtools'
include { MOSDEPTH                      } from '../modules/mosdepth'


workflow LONG_READ_VARIANT_CALLING {
	take:
		bam
		reference_genome
		reference_genome_index

	main:

        if (params.is_mapped) {
            MOSDEPTH(bam)
            coverage = MOSDEPTH.out
        } else {
            coverage = Channel.empty()
        }
		
		// Call SNPs
		if (params.call_snps) {
            // Split the BAM
		    split_bams = SPLITBAM(bam).split_bam.transpose()
			CLAIR3(split_bams, reference_genome, reference_genome_index)
			// snps = CLAIR3.out.gvcf
		} else {
			snps = Channel.empty()
		}

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
		// SNIFFLES2(bam, reference_genome, reference_genome_index)
		// SVIM(bam, reference_genome, reference_genome_index)
		// CUTESV(bam, reference_genome, reference_genome_index)
		// DYSGU(bam, reference_genome, reference_genome_index)
	
	emit:
		// sniffles  = SNIFFLES2.out.res_tuple
		// svim = SVIM.out.res_tuple
		// cutesv = CUTESV.out.res_tuple
		// dysgu = DYSGU.out.res_tuple
        coverage
        sniffles
        svim
        cutesv
        dysgu
}