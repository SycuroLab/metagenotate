#!/usr/bin/python
import os
import sys
import re
import csv
import argparse

parser = argparse.ArgumentParser()

# Sample command
# python prokka.py --genome_fasta_infile $genome_fasta_infile --num_cpus 40 --metagenome true --output_dir $sample_bin_prokka_annotation_dir;

fna_infile = None
metagenome = None
num_cpus = None
output_dir = None

parser.add_argument('--genome_fasta_infile', action='store', dest='fna_infile',
                    help='input file as input. (i.e. $HOME)')
parser.add_argument('--metagenome', action='store', dest='metagenome',
                    help='whether or not the genome is a MAG or not as input (i.e. true for MAGs)')
parser.add_argument('--num_cpus', action='store', dest='num_cpus',
                    help='number of cpus to run prokka. ')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='output directory as input. (i.e. $HOME)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

fna_infile = results.fna_infile
metagenome = results.metagenome
num_cpus = results.num_cpus
output_dir = results.output_dir

if(fna_infile == None):
	print('\n')
	print('error: please use the --genome_fasta_infile option to specify the input directory as input')
	print('fna_infile =' + ' ' + str(fna_infile))
	print('\n')
	parser.print_help()
	sys.exit(1)
if(metagenome == None):
	print('\n')
	print('error: please use the --metagenome option to specify the whether or not the genome is a MAG or not as input')
	print('metagenome =' + ' ' + str(metagenome))
	print('\n')
	parser.print_help()
	sys.exit(1)	
if(num_cpus == None):
	print('\n')
	print('error: please use the --num_cpus option to specify the number of cpus to run prokka as input')
	print('num_cpus =' + ' ' + str(num_cpus))
	print('\n')
	parser.print_help()
	sys.exit(1)
if(output_dir == None):
    print('\n')
    print('error: please use the --output_dir option to specify the output directory as input')
    print('output_dir =' + ' ' + str(output_dir))
    print('\n')
    parser.print_help()
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# The genome fasta file as input.
fasta_infile = fna_infile
	
# Get the fasta filename prefix.
fna_filename = os.path.basename(fna_infile)
fasta_filename = os.path.splitext(fna_filename)[0]

# If the input fna file is a gzipped file then uncompress the file.
if(fna_infile.endswith('.gz')):

	fasta_infile = os.path.splitext(fna_infile)[0]
	if(not(os.path.exists(fasta_infile))):

		os.system("gunzip -c {fna_infile} > {fasta_infile}".format(fna_infile=fna_infile,fasta_infile=fasta_infile))

annotations_dir = os.path.join(output_dir, "prokka_annotations")
if(not(os.path.exists(annotations_dir))):
	if(metagenome == "true"):
		os.system("prokka --metagenome --outdir {annotations_dir} --prefix {fasta_filename} {fasta_infile} --cpus {num_cpus} --rfam 1".format(annotations_dir=annotations_dir, fasta_filename=fasta_filename, fasta_infile=fasta_infile, num_cpus=num_cpus))
	else:
		os.system("prokka --outdir {annotations_dir} --prefix {fasta_filename} {fasta_infile} --cpus {num_cpus} --rfam 1".format(annotations_dir=annotations_dir, fasta_filename=fasta_filename, fasta_infile=fasta_infile, num_cpus=num_cpus))





'''
#### Going to use for gzipping files to save space.

# Tar compress annotations_dir if it exists and remove directory if tar compression was successful.
if(os.path.exists(annotations_dir)):
	os.system(("cd {output_dir}; tar -zcvf prokka_annotations.tar.gz prokka_annotations && rm -R prokka_annotations").format(output_dir=output_dir))

# Change permissions of output_dir to all read, write, and execute.
os.system("chmod -R 777 {output_dir}".format(output_dir=output_dir))

# Remove the uncompressed genome assembly genbank gbff file if it exists.
if(os.path.exists(genbank_infile)):
	os.system(("rm {genbank_infile}").format(genbank_infile=genbank_infile))
	
# Remove the uncompressed genome assembly fasta fna file if it exists.
if(os.path.exists(fasta_infile)):	
	os.system(("rm {fasta_infile}").format(fasta_infile=fasta_infile))

'''
