/*
 * Main Workflow: just_snps
 * 
 * Purpose:
 *   This is an archived standalone pipeline for SNP/small variant calling from long reads.
 *   It imports a dedicated workflow for long read SNP calling.
 *
 * Inputs (via parameters):
 *   - params.samplesheet: Sample sheet with sample IDs and BAM file paths
 *   - params.reference: Reference genome FASTA file
 *   - params.fai: Reference genome index file
 *   - params.minimap_index: Minimap2 index for the reference genome
 *
 * Other parameters:
 *   - params.model_path: Path to Clair3 model directory
 *   - params.sourceDir: Source directory for singularity container mounting
 *   - params.outdir: Output directory for results
 *
 * Process:
 *   1. Parse sample sheet to create input channel
 *   2. Call LONG_READ_VARIANTS workflow for variant calling
 *
 * Note:
 *   This is an archived pipeline that focuses only on SNP calling.
 */

nextflow.enable.dsl=2

/*
 * PARAMETERS
 */
params.samplesheet = null
params.reference = null
params.fai = null
params.model_path = "/hpcfs/users/a1767591/CONTAINERS/clair3_models/r1041_e82_400bps_sup_g615"
params.sourceDir = null
params.outdir = params.sourceDir


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
		Sample sheet:     ${params.samplesheet}
		Output directory: ${params.outdir}
		Reference genome: ${params.reference}
		Source directory: ${params.sourceDir}

=============================================================
"""

// Create input channel from the sample sheet
Channel.fromPath(params.samplesheet, checkIfExists: true)
					 .splitCsv(header: true)
		             .map ( row -> tuple(row.sampleID, row.bam, row.bai) )
					 .set { samples }


/*
 * Import the workflow from the workflows directory
 */
include { LONG_READ_VARIANTS as DLRVC } from './workflows/long_read_snps'

workflow {
	DLRVC(samples, params.reference, params.fai)
}


