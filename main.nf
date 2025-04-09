nextflow.enable.dsl=2

/*
 * PARAMETERS
 * Define all input parameters with default values where appropriate
 */
params.samplesheet = null  // Path to CSV containing sample information
params.reference = null    // Path to reference genome FASTA
params.reference_idx = null // Path to reference genome index (FAI)
params.minimap_index = null // Path to minimap2 index file (.mmi)
params.sourceDir = null    // Path to the desired results directory (needs to be already created prior to running the pipline. Needed by singularity)
params.outdir = params.sourceDir // Output directory (same as sourceDir by default)
params.skip_filtlong = false // Flag to skip filtlong filtering step
params.is_mapped = false   // Flag indicating if input reads are already mapped
params.call_snps = true    // Flag to enable SNP calling
params.sniffles = true     // Flag to enable Sniffles SV caller
params.svim = true         // Flag to enable SVIM SV caller
params.cutesv = true       // Flag to enable cuteSV SV caller
params.dysgu = true        // Flag to enable Dysgu SV caller

/*
 * CHECK INPUTS
 * Validate required input parameters and exit with error messages if missing
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

// Log the workflow parameters for user verification
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
        Pre-mapped data:         ${params.is_mapped}
        Call SNPs:               ${params.call_snps}
        Use Sniffles:            ${params.sniffles}
        Use SVIM:                ${params.svim}
        Use CuteSV:              ${params.cutesv}
        Use Dysgu:               ${params.dysgu}

=============================================================
"""

/*
 * INPUT CHANNEL SETUP
 * Create input channels from the sample sheet based on whether data is pre-mapped
 */
if (params.is_mapped) {
	log.info "Running variant calling only"
	// For pre-mapped data: Extract sampleID, fastq path (BAM), and BAM index
	Channel.fromPath(params.samplesheet, checkIfExists: true)
		   .splitCsv(header: true)
		   .map ( row -> tuple(row.sampleID, row.fastq, row.fastq + ".bai") )
		   .set { samples }
} else {
	log.info "Running preprocessing, mapping and variant calling"
	// For unmapped data: Extract sampleID and bam or fastq path
	Channel.fromPath(params.samplesheet, checkIfExists: true)
		   .splitCsv(header: true)
		   .map ( row -> tuple(row.sampleID, row.fastq) )
		   .set { samples }
}


/*
 * IMPORT WORKFLOWS
 * Include workflow modules from external files
 */
include { LONG_READ_VARIANTS as DLRVC      } from './workflows/long_read_variants'  // Main workflow for unmapped reads
include { SAMTOOLS_FAIDX                   } from './modules/samtools'  // Samtools faidx module
// include { LONG_READ_SV_CALLING as ONLY_SVS } from './workflows/long_read_svs_only'  // SV calling workflow (commented out)
include { ALREADY_MAPPED as VARIANTS       } from './workflows/already_mapped'  // Workflow for pre-mapped reads
// include { LONG_READ_VARIANTS as SNPS 	   } from './workflows/long_read_snps'  // SNP calling workflow (commented out)

/*
 * WORKFLOW EXECUTION
 * Define the entry point based on input read status
 */
workflow {
    if (params.is_mapped) {
        // If the data is already mapped, run the VARIANTS workflow
        VARIANTS(samples)
    }
    else {
        // If the data is not mapped, run the full DLRVC workflow
		DLRVC(samples)
	}
}


