# DaviesCentreInformatics/LR-variantCaller: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) 
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- List of new features

### Changed

- List of what was changed

### Deprecated

- List of what will soon be removed

### Removed

- List of what no longer works

### Fixed

- List of bugs that were fixed

### Security

- List of security concerns/vulnerabilities

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
