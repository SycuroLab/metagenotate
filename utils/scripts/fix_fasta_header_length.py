#!/usr/bin/python
import os
import sys
import re
import csv
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()

fasta_infile = None
fasta_outfile = None

parser.add_argument('-i', action='store', dest='fasta_infile',
                    help='input fasta file as input. (i.e. infile.fasta)')
parser.add_argument('-o', action='store', dest='fasta_outfile',
                    help='output fasta file as output. (i.e. outfile.fasta)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

fasta_infile = results.fasta_infile
fasta_outfile = results.fasta_outfile

if(fasta_infile == None):
	print('\n')
	print('error: please use the -i option to specify the input fasta file as input')
	print('fasta_infile =' + ' ' + str(fasta_infile))
	print('\n')
	parser.print_help()
	sys.exit(1)
if(fasta_outfile == None):
    print('\n')
    print('error: please use the -o option to specify the output fasta file as output')
    print('fasta_outfile =' + ' ' + str(fasta_outfile))
    print('\n')
    parser.print_help()
    sys.exit(1)

basename = os.path.basename(fasta_outfile)
filename = os.path.splitext(basename)[0]

# Get the output directory from the output fasta file and make the directory if it doesn't exist.
output_dir = os.path.dirname(fasta_outfile)
if not os.path.exists(output_dir):
	os.makedirs(output_dir)


fasta_output_file = open(fasta_outfile, "w+")
for record in SeqIO.parse(fasta_infile, "fasta"):
    seq_id = record.id
    desc = record.description
    seq = record.seq
	
    print(len(seq_id))
    
    # Fix the header so that the header length is less than 37 characters in length for prokka.
    # The rest of the header for metaspades is usually decimal places in the coverage calculation.
    if(len(seq_id) >= 37):
        print(seq_id)
        print(desc)
        seq_id = seq_id[0:36]
        seq_length = len(seq)
        
    # Make a new record and print to a new file.
    record = SeqRecord(
        Seq(str(seq)),
        id=seq_id,
        name="",
        description=desc,
    )
    SeqIO.write(record, fasta_output_file, "fasta")
fasta_output_file.close()

