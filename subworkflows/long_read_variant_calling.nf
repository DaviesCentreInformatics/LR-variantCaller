/*
 * Workflow: LONG_READ_VARIANT_CALLING
 * 
 * Purpose:
 *   This workflow performs variant calling on long read data, including both SNPs/small variants 
 *   and structural variants (SVs) using various tools.
 *
 * Inputs:
 *   - bam: Channel containing BAM files with long reads aligned to reference
 *   - reference_genome: Reference genome FASTA file
 *   - reference_genome_index: Index file for the reference genome
 *
 * Outputs:
 *   - coverage: Sequencing coverage information from MOSDEPTH
 *   - sniffles: Structural variants called by Sniffles2 (optional)
 *   - svim: Structural variants called by SVIM (optional)
 *   - cutesv: Structural variants called by cuteSV (optional)
 *   - dysgu: Structural variants called by Dysgu (optional)
 *
 * Parameters:
 *   - params.is_mapped: Flag indicating if reads are already mapped
 *   - params.call_snps: Flag to enable SNP/small variant calling
 *   - params.sniffles: Flag to enable Sniffles2 SV calling
 *   - params.svim: Flag to enable SVIM SV calling
 *   - params.cutesv: Flag to enable cuteSV SV calling
 *   - params.dysgu: Flag to enable Dysgu SV calling
 */

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
        // Calculate coverage statistics if reads are already mapped
        if (params.is_mapped) {
            MOSDEPTH(bam)
            coverage = MOSDEPTH.out
        } else {
            coverage = Channel.empty()
        }
		
		// Call SNPs using Clair3 when enabled
		if (params.call_snps) {
            // Split the BAM file by chromosome for parallel processing
		    split_bams = SPLITBAM(bam).split_bam.transpose()
			CLAIR3(split_bams, reference_genome, reference_genome_index)
			// snps = CLAIR3.out.gvcf
		} else {
			snps = Channel.empty()
		}

		// Call structural variants using multiple tools based on parameters
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