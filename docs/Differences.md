# Overview of nf-EXPLOR and ours

## 1. QC

Both use `filtlong` to perform QC. Both use `NanoPlot` to generate before and after QC plots.

`filtlong` will only work with fastq files so first we convert the bam to fastq then run `filtlong` on the fastq files then take the readnames of reads that pass the filtering and extract those read names from the bam file.

## 2. Mapping

`nf-EXPLOR` uses minimap2 to map reads to the reference genome. We use `dorado` to map modbams to the reference. `dorado` used minimap2 to perform the alignments

## 3. SNP calling

Both use `Clair3`. The main difference being the model used to call SNPs. Ours matches the model we used to basecall the reads.

## 4. SV calling

We use `sniffles2`, `dysgu`, `svim` and `cuteSV` to call SVs. `nf-EXPLOR` uses `sniffles` only, by default. Where possible, I've used the same parameters as `nf-EXPLOR`.