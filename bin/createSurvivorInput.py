#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import sys

# Parse command line arguments
parser = argparse.ArgumentParser(description='Create input files for Survivor.')
parser.add_argument('--sampleID', type=str, default='sampleID',
					help='Name/ID of the sample')
parser.add_argument('--output', type=str, default='survivor_input.txt',
					help='Name of the file that will contain the file paths of '
					'input files.')
parser.add_argument('vcfs', type=str, nargs='+', help='vcf files to be processed')
args = parser.parse_args()

# Unzip the vcf files
for vcf in args.vcfs:
	if vcf.endswith('.gz'):
		os.system(f'unpigz -f {vcf}')

vcfs = [os.path.splitext(vcf)[0] for vcf in args.vcfs]

df = pd.DataFrame({'Path': vcfs})
df['Sample'] = args.sampleID
df['Chrom'] = df['Path'].apply(lambda x: os.path.basename(x)
							   .split('.')[0]
							   .split('_')[-1])
df.sort_values(by='Chrom', inplace=True)
outpath = os.path.join(os.getcwd(), f'{args.sampleID}_{args.output}')
df['Path'].to_csv(outpath, index=False, header=False, sep='\t')