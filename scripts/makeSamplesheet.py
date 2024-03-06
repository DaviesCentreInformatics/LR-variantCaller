#!/usr/bin/env python3

import os
import re
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(
	description='Create a sample sheet for the SNP calling\
pipeline\n\nExample:\npython makeSamplesheet.py -i \
/path/to/input -o /path/to/output/samplesheet.csv',
	formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument('-i', '--input',
					help='Complete path to the input directory',
					required=True)

parser.add_argument('-o', '--output',
					help='Path to the output samplesheet',
					required=True)

args = parser.parse_args()

#print(f"Input directory: {args.input}")
#print(f"Output samplesheet: {args.output}")


# Check parameters provided are valid
if not os.path.abspath(args.input):
	print(f"Error: {args.input} is not an absolute path")
	exit(1)

if not os.path.exists(os.path.dirname(args.output)):
	print(f"Error: {args.output} does not exist")
	exit(1)

if os.path.exists(args.output):
	print(f"Error: {args.output} already exists")
	exit(1)

# List the contents of the input path
input_files = os.listdir(args.input)

# Filter the list of files to only include fastq files
fastq_files = [f for f in input_files if f.endswith(
	('fq', 'fastq', 'fastq.gz', 'fq.gz'))]

# Check that `fastq_files` is not empty
if len(fastq_files) == 0:
	print(f"Error: No fastq files found in {args.input}")
	print("Perhaps they end in something other than \
.fq, .fastq, .fq.gz, or .fastq.gz")
	exit(1)

# Filter so that we have two lists, one for the forward reads and one for the reverse reads
forward_reads = sorted(
	[f for f in fastq_files if re.search(r'.*(_1|_R1).*', f)]
	)
reverse_reads = sorted(
	[f for f in fastq_files if re.search(r'.*(_2|_R2).*', f)]
	)

# Generate the sample names
sample_names = [f.split('_')[0] for f in forward_reads]

assert len(forward_reads) == len(reverse_reads), "Error: Unequal number of forward and reverse reads"
assert len(sample_names) == len(set(sample_names)), "Error: Duplicate sample names detected"
assert len(sample_names) == len(forward_reads), "Error: Unequal number of sample names and forward reads"

input_path = args.input

print('Writing sample sheet to {}'.format(args.output))
with open(args.output, 'w') as f:
	f.write('sampleID,read1,read2\n')
	for sample, read1, read2 in zip(sample_names, forward_reads, reverse_reads):
		f.write(f'{sample},{os.path.join(input_path,read1)},{os.path.join(input_path,read2)}\n')

print(f"Sample sheet written to {args.output}")
