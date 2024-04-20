# DaviesCentreInformatics/LR-variantCaller: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) 
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `coverage.nf` as a "subworkflow" to just run `NanoPlot` on raw ONT reads.
- Added `--skip_filtlong` flag to skip read filtering. 
  Default is to run `FiltLong`.
- 

### Changed

- `clair3.nf` now has output specified. This should now report instances where 
  the output is not generated.
- Reformatted `README.md` to be follow a template and be more informative.
- `clair3` now uses v1.0.7 and has new ONT models added.

### Deprecated

- N/A

### Removed

- N/A

### Fixed

- N/A

### Security

- N/A

## [v1.0.2](https://github.com/DaviesCentreInformatics/LR-variantCaller/releases/tag/v1.0.2)

### Fixed

- Fixed error where splitbam would split all contigs, not just actual chromosomes.

## [v1.0.1](https://github.com/DaviesCentreInformatics/LR-variantCaller/releases/tag/v1.0.1) - 2024-03-07

### Added

- `makeSamplesheet.py` script to generate a samplesheet from a directory of fastq files.

### Fixed

- Fixed mismatch in output pattern and input patterns for `sniffles2.nf` after changing the `splitbam` output pattern.

## [v1.0.0](https://github.com/DaviesCentreInformatics/LR-variantCaller/releases/tag/v1.0.0) - 2024-03-06

### Initial release

The initial release of the Davies Centre Informatics Oxford Nanopore Technology Structural Variant Calling Pipeline.
The pipeline performs the following steps:

1. QC of raw reads with NanoPlot
2. Long read trimming with FiltLong
3. QC of trimmed reads with NanoPlot
4. Mapping to reference genome with Minimap2
5. Analyse depth and coverage with MosDepth
6. Flag, Index and general alignment stats with SAMtools
7. Split bams per chromosome per sample.
8. Identify SNP and SV variants per sample per chromosome
9. Generate MultiQC report.

**Full Changelog**: https://github.com/DaviesCentreInformatics/LR-variantCaller/commits/v1.0.0
 