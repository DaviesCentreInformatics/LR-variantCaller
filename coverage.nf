nextflow.enable.dsl=2

/*
 * PARAMETERS
 */
params.samplesheet = null
params.sourceDir = null
params.outdir = params.sourceDir


/*
 * CHECK INPUTS
 */
if (params.samplesheet == null) {
	error "Please provide a sample sheet using the `--samplesheet` flag"
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
N E X T F L O W  --  C O V E R A G E -- O N T
=============================================================

Pipeline Parameters:
		Sample sheet:     ${params.samplesheet}
		Output directory: ${params.outdir}
		Source directory: ${params.sourceDir}

=============================================================
"""

// Create input channel from the sample sheet
Channel.fromPath(params.samplesheet, checkIfExists: true)
					 .splitCsv(header: true)
		             .map ( row -> tuple(row.sampleID, row.fastq) )
					 .set { samples }

include { NANOPLOT } from './modules/nanoplot'

workflow {
	/*
	 * NANOplot
	 */
	NANOPLOT(samples)
}