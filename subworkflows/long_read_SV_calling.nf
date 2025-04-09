/*
 * Workflow: LONG_READ_SV_CALLING
 * 
 * Purpose:
 *   This workflow focuses specifically on calling structural variants (SVs) from long read data
 *   using multiple SV callers to enable ensemble or comparative approaches.
 *
 * Inputs:
 *   - mapped: Channel containing aligned BAM files
 *   - reference_genome: Reference genome FASTA file
 *   - reference_genome_index: Index file for the reference genome
 *
 * Outputs:
 *   - sniffles: Structural variants called by Sniffles2 (optional)
 *   - svim: Structural variants called by SVIM (optional)
 *   - cutesv: Structural variants called by cuteSV (optional)
 *   - dysgu: Structural variants called by Dysgu (optional)
 *
 * Parameters:
 *   - params.sniffles: Flag to enable Sniffles2 SV calling
 *   - params.svim: Flag to enable SVIM SV calling
 *   - params.cutesv: Flag to enable cuteSV SV calling
 *   - params.dysgu: Flag to enable Dysgu SV calling
 */

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
		// Call SVs with different tools based on parameter flags
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