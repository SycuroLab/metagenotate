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
                    help='input fasta file path as input. (i.e. sequences.fasta)')
parser.add_argument('-o', action='store', dest='fasta_outfile',
                    help='output fasta file path as input. (i.e. sequences.fasta)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

fasta_infile = results.fasta_infile
fasta_outfile = results.fasta_outfile

if(fasta_infile == None):
	print('\n')
	print('error: please use the -i option to specify the input fasta file path as input')
	print('fasta_infile =' + ' ' + str(fasta_infile))
	print('\n')
	parser.print_help()
	sys.exit(1)
if(fasta_outfile == None):
    print('\n')
    print('error: please use the -o option to specify the output fasta file path as input')
    print('fasta_outfile =' + ' ' + str(fasta_outfile))
    print('\n')
    parser.print_help()
    sys.exit(1)


basename = os.path.basename(fasta_outfile)
filename = os.path.splitext(basename)[0]
output_dir = os.path.dirname(fasta_outfile)

if not os.path.exists(output_dir):
	os.makedirs(output_dir)

sequence_length_dict = {}
for record in SeqIO.parse(fasta_infile, "fasta"):
    seq_id = record.id
    desc = record.description
    seq = record.seq
    seq_length = len(seq)
    
    if(not(seq_length in sequence_length_dict)):
        sequence_length_dict[seq_length] = []
    if(seq_length in sequence_length_dict):
        sequence_length_dict[seq_length].append(record)

fasta_output_file = open(fasta_outfile, "w+")
for seq_length in sorted(sequence_length_dict):
    for record in sequence_length_dict[seq_length]:
        SeqIO.write(record, fasta_output_file, "fasta")
fasta_output_file.close()



