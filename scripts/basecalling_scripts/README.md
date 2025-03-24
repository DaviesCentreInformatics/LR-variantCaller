# Base calling POD5 on Phoenix

Here are the general steps for base calling POD5 files to UBAM files.

1. Convert POD5 files to UBAMs with `1_pod5_to_ubam.sh`. I found the best way
   to do this was to give it the relatively small individual POD5 files and then
   let it convert each one rather than merging into a big POD5 file. This way
   the conversion happens quite quickly and you can get jobs submitted quickly
   because SLURM sees you're only requesting the GPU for an hour or two.

``` bash
for pod5 in /path/to/pod5/files/*.pod5;
do
    # sbatch 1_pod5_to_ubam.sh pod5_file output_dir model
    sbatch 1_pod5_to_ubam.sh $pod5 $output_dir sup
done
```

2. Check the integrity of each UBAM before concatenating them together.
   
``` bash
samtools quickcheck -u /path/to/ubam/files/*.bam
```

Use `-u` to allow `samtools` to ignore the fact that these are unmapped BAM
files.

3. Concatenate the UBAM files together with `samtools cat`.

``` bash
# Create a list of input files.
find /path/to/ubams -name "*.bam" > ubam_files.txt

# Concatenate the files together.
samtools cat -@ 16 -b ubam_files.txt -o /path/to/ubams/concatenated.bam

# Check the integrity before deleting the individual UBAM files.
samtools quickcheck -u /path/to/ubams/concatenated.bam
```

4. Demultiplex the concatenated UBAM file with `2_demux_ubam.sh`
