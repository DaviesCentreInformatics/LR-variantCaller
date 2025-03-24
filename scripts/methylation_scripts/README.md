# Steps for calling methylation values with Modkit

1. Map the UBAM files to the reference genome with minimap2. You **MUST** ensure
   the methyl tags are presevered in the bam file prior to mapping or it will
   fail. Use `1_DORADO_ALIGN.sh` to do this.
2. Sort the BAM files with `samtools sort`. Use `2_DORADO_SORT.sh` to do this.
3. Call the methylation values with `modkit pileup`. Use `3_MODKIT_PILEUP.sh` to
   do this.

