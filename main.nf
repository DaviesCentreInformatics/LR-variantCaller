nextflow.enable.dsl=2

/*
 * PARAMETERS
 */
params.samplesheet = null
params.reference = null
params.reference_idx = null
params.minimap_index = null
params.sourceDir = null
params.outdir = params.sourceDir
params.skip_filtlong = false
params.already_mapped = null


/*
 * CHECK INPUTS
 */

if (params.samplesheet == null) {
	error "Please provide a sample sheet using the `--samplesheet` flag"
	System.exit(1)
}
if (params.minimap_index == null) {
	error "Please provide a minimap_index using the `--minimap_index` flag"
	System.exit(1)
}

if (params.reference == null) {
	error "Please provide a reference genome using the `--reference` flag"
	System.exit(1)
}
if (params.reference_idx == null) {
	error "Please provide a reference genome index using the `--reference_idx` flag"
	System.exit(1)
}

if (params.sourceDir == null) {
	error """Please provide a source directory using the `--sourceDir` flag
	         This should be the identical to the outdir and is used by the
			 singularity container to mount the results directory"""
	System.exit(1)
}

if (params.outdir != params.sourceDir) {
	error "The outdir and sourceDir parameters must be identical"
	System.exit(1)
}

log.info """\
N E X T F L O W  --  V A R I A N T   C A L L I N G  --  O N T
=============================================================

Pipeline Parameters:
		Sample sheet:            ${params.samplesheet}
		Reference genome:        ${params.reference}
		Reference genome index:  ${params.reference_idx}
		Minimap2 index:          ${params.minimap_index}
		Clair3 model:            ${params.model_path} 
		Source/Output directory: ${params.sourceDir}
		Skip SNP calling:        ${params.skip_snps}
		Methylation:             ${params.methylation}
		Only SVs:                ${params.only_svs}

=============================================================
"""

// Create input channel from the sample sheet
if (params.only_svs || params.already_mapped) {
	log.info "Running SV calling only"
	Channel.fromPath(params.samplesheet, checkIfExists: true)
					 .splitCsv(header: true)
		             .map ( row -> tuple(row.sampleID, row.fastq, row.fastq + ".bai") )
					 .set { samples }
} else {
	log.info "Running variant calling"
	Channel.fromPath(params.samplesheet, checkIfExists: true)
					 .splitCsv(header: true)
		             .map ( row -> tuple(row.sampleID, row.fastq) )
					 .set { samples }
}


/*
 * Import the workflow from the workflows directory
 */
include { LONG_READ_VARIANTS as DLRVC      } from './workflows/long_read_variants'
include { SAMTOOLS_FAIDX                   } from './modules/samtools'
include { LONG_READ_SV_CALLING as ONLY_SVS } from './workflows/long_read_svs_only'
include { ALREADY_MAPPED as VARIANTS       } from './workflows/already_mapped'

workflow {
	if (params.only_svs) {
		
		// (fasta, fai) = SAMTOOLS_FAIDX(params.reference)
		fasta = params.reference
		fai = params.reference_idx

		params.sniffles = true
		params.svim = false
		params.cutesv = true
		params.dysgu = false
		ONLY_SVS(samples, fasta, fai)
	
	} else if (params.already_mapped) {
		VARIANTS(samples)
		
	} else {
		DLRVC(samples)
	}
}


