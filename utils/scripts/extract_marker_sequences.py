#!/usr/bin/python
import os
import sys
import re
import csv
import argparse
import glob

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# conda activate extract_marker_sequences_env
# python extract_marker_sequences.py --fasta_infile ~/Desktop/macbook_air/graduate_project/2021_pullulanase_paper/prokka_annotation_dir/prokka_annotation_dir/GCF_009730275.1_ASM973027v1_genomic.fna.fna --gff_infile ~/Desktop/macbook_air/graduate_project/2021_pullulanase_paper/prokka_annotation_dir/prokka_annotation_dir/GCF_009730275.1_ASM973027v1_genomic.fna.gff --output_dir ~/Desktop/macbook_air/graduate_project/2021_pullulanase_paper/prokka_annotation_dir

parser = argparse.ArgumentParser()

fasta_infile = None
gff_infile = None
output_dir = None

parser.add_argument('--fasta_infile', action='store', dest='fasta_infile',
                    help='The genome assembly file, metagenome assembly file or metagenome bin file as input. (i.e. filename.fasta)')
parser.add_argument('--gff_infile', action='store', dest='gff_infile',
                    help='The gff3 formated file generated from an annotation program as input. (i.e. filename.gff)')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='The output directory as input. (i.e. $HOME)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

fasta_infile = results.fasta_infile
gff_infile = results.gff_infile
output_dir = results.output_dir

if(fasta_infile == None):
	print('\n')
	print('error: please use the --fasta_infile option to specify the input directory as input')
	print('fasta_infile =' + ' ' + str(fasta_infile))
	print('\n')
	parser.print_help()
	sys.exit(1)
if(gff_infile == None):
	print('\n')
	print('error: please use the --gff_infile option to specify the input directory as input')
	print('gff_infile =' + ' ' + str(gff_infile))
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


# Marker information. Marker name, search term text, and sequence type.
markers = { 
	'cpn60': {
			'search_terms': ['60 kDa chaperonin'],
			'type': 'CDS'
	},
	'16S': {
			'search_terms': ['16S ribosomal RNA'],
			'type': 'rRNA'
	}	
}


fna_filename = os.path.basename(fasta_infile)

fasta_filename = os.path.splitext(fna_filename)[0]

fasta_file = os.path.splitext(fasta_infile)[0]

print(fasta_file)

# Open genome assembly fasta file as a file handle for reading.
fasta_input_file = open(fasta_infile)

# Parse genome assembly fasta file and construct a sequence dictionary where sequence ids as key and sequence as the value.
sequence_dict = SeqIO.to_dict(SeqIO.parse(fasta_input_file, "fasta"))

# Close the fasta input fasta file handle.
fasta_input_file.close()

# Open the GFF annotation file as a file handle for reading.
gff_input_file = open(gff_infile, "r+")

# Dictionary for sequence entries parsed from the genome assembly file and GFF file.
sequence_entries = {}

## GFF input file fields.

#Fields
#Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'
#
#seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
#source - name of the program that generated this feature, or the data source (database or project name)
#feature - feature type name, e.g. Gene, Variation, Similarity
#start - Start position* of the feature, with sequence numbering starting at 1.
#end - End position* of the feature, with sequence numbering starting at 1.
#score - A floating point value.
#strand - defined as + (forward) or - (reverse).
#frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

for line in gff_input_file.readlines():
    #print(line)
    
    # Break out of loop if "##FASTA".
    if(("##FASTA" in line)):
        break
        
    # If not comment lines and GFF file content.
    if(not("##" in line)):
        
        line = line.strip("\n")
        (sequence_id, source, feature_type, start, end, score, strand, frame, attribute) = line.split("\t")
        #NZ_CP046311.1 Prodigal:002006 CDS 533 2956 . - 0 ID=FLKKEMPL_00001;Name=spoIIIE;db_xref=COG:COG1674;gene=spoIIIE;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P21458;feature_id=FLKKEMPL_00001;product=DNA translocase SpoIIIE
        #print(sequence_id, source, feature_type, start, end, score, strand, frame, attribute)
        
        #print(feature_type)
        #print(attribute)
        
        #ID=FLKKEMPL_00001;Name=spoIIIE;db_xref=COG:COG1674;gene=spoIIIE;inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P21458;feature_id=FLKKEMPL_00001;product=DNA translocase SpoIIIE
        # Split by ';' delimiter and make a dictionary based on the attribute names.
        attributes = attribute.split(";")
        attributes_metadata = {}
        for attribute in attributes:
            # split attributes by the '=' symbol delimiter.
            (key,value) = attribute.split("=")
            attributes_metadata[str(key)] = value
        #print(attributes_metadata)
        
        
        
        for marker in markers:

            search_terms = markers[marker]["search_terms"]
            type = markers[marker]["type"]

            if(feature_type == type):

                for search_term in search_terms:
                    #print(search_term)
                    if(search_term in str(attributes_metadata["product"])):

                        feature_metadata = attributes_metadata
                        feature_id = feature_metadata["ID"]
                        
                        feature_metadata["type"] = feature_type
                        feature_metadata["location"] = {"start": int(start), "end": int(end), "strand": strand}
                        #print(feature_metadata)
                        #print(start,end)

                        # Initialize the start and end index values.
                        start_index = 0
                        end_index = 0

                        if(strand == "+"):

                            # Get the start position index.
                            start_index = (int(start) - 1)

                            # Get the end position index.
                            end_index = int(end)

                        elif(strand == "-"):

                            # Get the start position index.
                            start_index = (int(start) - 1)

                            # Get the end position index.
                            end_index = (int(end))

                        # Get the sequence string by sequence id from the sequence dictionary.
                        feature_sequence = str(sequence_dict[sequence_id][int(start_index):int(end_index)].seq)
                        #feature_sequence = sequence[int(start):int(end)]
#                        print(strand)
                        if(strand == "-"):
                            feature_sequence = str(Seq(feature_sequence).reverse_complement())
                        #print(feature_sequence)
                        #print("{}".format(sequence_id))
                        sequence_length = len(feature_sequence)
                        #print(feature_sequence)
                        #print(sequence_length)
                        
                        if(not(marker in sequence_entries)):
                            sequence_entries[marker] = {}
                        if(not(sequence_id in sequence_entries[marker])):
                            sequence_entries[marker][sequence_id] = []
                            sequence_entries[marker][sequence_id].append([feature_id,feature_metadata,search_term,feature_type,sequence_length,feature_sequence])
                        else:
                            sequence_entries[marker][sequence_id].append([feature_id,feature_metadata,search_term,type,sequence_length,feature_sequence])

# Print out marker sequence metadata file for parsing.
for marker in markers:
    csv_writer_file_handle = open(os.path.join(output_dir, "_".join([marker,"metadata"]) + ".csv"), "w+")
    csv_writer = csv.writer(csv_writer_file_handle, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
    csv_writer.writerow(["sequence_id","feature_id","marker","search_term","type","feature_metadata","marker_sequence"])
    if(marker in sequence_entries):
        fasta_outfile = os.path.join(output_dir, "_".join([marker, "sequences"]) + ".fasta")
        fasta_output_file = open(fasta_outfile, "w+")
        for sequence_id in sequence_entries[marker]:
            for sequence_entry in sequence_entries[marker][sequence_id]:
                (feature_id,feature_metadata,search_term,type,sequence_length,marker_sequence) = sequence_entry
                
                csv_writer.writerow([sequence_id,feature_id,marker,search_term,type,feature_metadata,marker_sequence])
                
                source = "Prokka"
                sequence_length = len(marker_sequence)

                header = "{feature_id} scaffold_id={sequence_id} locus_tag={feature_id} marker={marker} source={source} length={sequence_length}".format(sequence_id=sequence_id,feature_id=feature_id,sequence_length=sequence_length,marker=marker,source=source)
#                print(header)
                record = SeqRecord(
                    Seq(str(marker_sequence)),
                    id=header,
                    description=""
                )


                SeqIO.write(record, fasta_output_file, "fasta")
        fasta_output_file.close()


