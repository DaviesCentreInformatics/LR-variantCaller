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

3. Check the methyl tags are present in the UBAM files.

``` bash
samtools view -h /path/to/ubam/files/a_single_ubam.bam | less
```

``` console
@HD     VN:1.6  SO:unknown
@PG     ID:basecaller   PN:dorado       VN:0.5.3+d9af343        CL:dorado basecaller -v -b 0 --min-qscore 9 /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/SOFTWARE/dorado-0.5.3-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 --modified-bases-models /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/SOFTWARE/dorado-0.5.3-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mCG_5hmCG@v2/ /hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/TEST_DATA/ONT/POD5/20240326_1335_P2S-01629-A_PAW21090_76f4a9ea/pod5_pass/barcode01/PAW21090_pass_barcode01_76f4a9ea_a02d5188_0.pod5
@PG     ID:samtools     PN:samtools     PP:basecaller   VN:1.17 CL:samtools view -h PAW21090_pass_barcode01_76f4a9ea_a02d5188_0_sup.bam
@RG     ID:a02d5188d93a9db31db03aa20c8a12a50e9a784d_dna_r10.4.1_e8.2_400bps_sup@v4.2.0  PU:PAW21090     PM:scylla       DT:2024-03-26T03:07:42.401+00:00        PL:ONT  DS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v4.2.0 modbase_models=dna_r10.4.1_e8.2_400bps_sup@v4.2.0_5mCG_5hmCG@v2 runid=a02d5188d93a9db31db03aa20c8a12a50e9a784d     LB:Sample1-3    SM:Sample1-3
8a484171-0733-4d78-bc54-564926c5ced9    4       *       0       0       *       *       0       0       TATGTTAACCTACTGGTTCAGTTGCATGTGTAAGGTTAACACAAAGACACCAACAACTTTCTCAGCACCTTCCTGCTTTTGATTACCCTGACTTTTCATTGTATTTTATAATCTTCATTGCTTTTATTACTCAGTACTTCCCATCTGTTTTCATAATTTTTTCCCCAGGAAAGATTCAAGATCATTTTTATCCAAGGGTGAAGTTGAGAATAATTTAGAAGACTTAACATTTGACAATGATAAGAATTTTTCTTATTCAGCTCTTCAGGAAAATATTATAGACTCTTGCACA    $%''('$$%%$$$%&%')(),10('%%%&&(''(245546660000((&&%%%%')('*('('()56655679:<=>>=868+*,10.-***))*.(/79877778:9765567888:::::;;;1..)&&%%&('&&'&&+'(69999:87886540++++*-..09998966666886656689;=>:744447842''''&&%%'())(,,,,../01298522399:::::9887666667898;>>=;96668777788:)())''4-.+**,--+++++,,//'''    qs:i:12 du:f:0.7264     ns:i:3632       ts:i:10 mx:i:1  ch:i:354        st:Z:2024-03-26T03:08:48.253+00:00      rn:i:67 fn:Z:PAW21090_pass_barcode01_76f4a9ea_a02d5188_0.pod5   sm:f:102.174    sd:f:27.7612    sv:Z:quantile   dx:i:0  RG:Z:a02d5188d93a9db31db03aa20c8a12a50e9a784d_dna_r10.4.1_e8.2_400bps_sup@v4.2.0        MN:i:292
        MM:Z:C+h?;C+m?; ML:B:C
```

The `@PG` line shows the command that was used to generate the BAM file.  
The `SM:` tag in the `@RG` line shows the sample name. As this was a multiplexed
sample, the sample name is `Sample1-3`.  
The `MM:` and `ML:` tags describe the methylation information. For details on
how to interpret these tags, see pages 7-9 from [here](https://samtools.github.io/hts-specs/SAMtags.pdf).

4. Concatenate the UBAM files together with `samtools cat`.

``` bash
# Create a list of input files.
find /path/to/ubams -name "*.bam" > ubam_files.txt

# Concatenate the files together.
samtools cat -@ 16 -b ubam_files.txt -o /path/to/ubams/concatenated.bam

# Check the integrity before deleting the individual UBAM files.
samtools quickcheck -u /path/to/ubams/concatenated.bam
```

5. Demultiplex the concatenated UBAM file with `2_demux_ubam.sh`
