#!/usr/bin/env python3
###################################
# File      - DNAcleaner
# Modified  - sat 28 maj 2022 15:37:18 CEST
# Sign      - Abhijeet
###################################
program = 'DNAcleaner'
version = ': version (1.0)\nAuthor: Abhijeet Singh <abhijeetsingh.aau@gmail.com>'
citation = '''\n\nCitation: Singh, Abhijeet. DNAcleaner: a python utility to clean DNA sequences and headers.
ResearchGate 2022, http://dx.doi.org/10.13140/RG.2.2.30762.29124/1, Available at GitHub: https://github.com/abhijeetsingh1704/DNAcleaner'''
###################################
import sys
import os
import datetime
import subprocess
from collections import defaultdict
import argparse
import string
import re
#
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Biopython missing, attempting to install...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython>=1.79"])
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
#
###################################
# parser
parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter, 
    description='Cleans invalid bases/residues in DNA sequence file i.e., replaces invalid bases with Ns and optionally removes special characters from headers' + citation)
optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--input", dest='input', required=True, type=str, help='input file')
optional.add_argument("-f", "--input_format", dest='input_format', required=False, type=int, default=1, help='format of input file\n1 = fasta (default)\n2 = genbank\n')
optional.add_argument("-o", "--output", dest='output', required=False, type=str, help='output file (default: Clean<time>_<input_file><.ext>)')
optional.add_argument("-F", "--output_format", dest='output_format', required=False, type=int, default=1, help='format of output fasta file\n1 = interleaved fasta (default)\n2 = fasta-2line\n3 = tab-delimited\n')
optional.add_argument("-c", "--clean_header", dest='clean_header', required=False, metavar='Y/y or N/n', default='N', type=str, help='clean and replace special characters in header\nwith underscore "_" (default: N)')
optional.add_argument("-v", "--verbose", dest='verbose', required=False, metavar='Y/y or N/n', default='Y', type=str, help='print progress to the terminal (default: verbose)')
optional.add_argument('-V', '--version', action='version', version='%(prog)s'+ str(version))
parser._action_groups.append(optional)
args = parser.parse_args()
###################################
#
def filenames(filename):
    absname = os.path.abspath(filename)
    basename = os.path.basename(filename)
    file_name = os.path.splitext(basename)[0]
    file_ext = os.path.splitext(basename)[1]
    return absname, basename, file_name, file_ext
###################################
#
date_var_1 = datetime.datetime.now()
date_var_2 = date_var_1.strftime("%H%M%S")
###################################
# 
input_args = args.input
output_args = args.output
input_format_args = args.input_format
output_format_args = args.output_format
verbosity_args = args.verbose.upper()
cleaner_args = args.clean_header.upper()
#
if cleaner_args == "Y":
    clean_word = "YES"
elif cleaner_args == "N":
    clean_word = "NO"
###################################
#
input_file = filenames(input_args)[0]
###################################
# optional options input format
if type(input_format_args) == int:
    if input_format_args == 1:
        input_format = "fasta"
    if input_format_args == 2:
        input_format = "genbank"
    if input_format_args >= 3:
        sys.exit("[ERROR]: input format flag requires interger value\n\t use -f or --input_format with value 1 or 2")
else:
    sys.exit("[ERROR]: input format flag requires interger value\n\t use -f or --input_format with value 1 or 2")
###################################
# optional options output format
if type(output_format_args) == int:
    if output_format_args == 1:
        output_format = "fasta"
        output_ext = ".fasta"
    if output_format_args == 2:
        output_format = "fasta-2line"
        output_ext = ".fasta"
    if output_format_args == 3:
        output_format = "tab"
        output_ext = ".tab"
    if output_format_args > 3:
        sys.exit("[ERROR]: output format flag requires interger value\n\t use -F or --output_format with value 1, 2 or 3")
else:
    sys.exit("[ERROR]: output format flag requires interger value\n\t use -F or --output_format with value 1, 2 or 3")
###################################
# output
if output_args is not None:
    output_var = os.path.abspath(output_args)
    output_file = output_var + output_ext
    output_obj = open(output_file, "w")
else:
    output_var = os.path.join(os.getcwd(), "Clean" + date_var_2 + "_" + filenames(input_file)[2])
    output_file = output_var + output_ext
    output_obj = open(output_file, "w")
###################################
# change objects
DNA_bases = list("ATGCN")
numbers = str("0123456789")
All_Characters = list(string.ascii_uppercase + string.punctuation + numbers) #.replace(".","")
Remove_Characters = [x for x in All_Characters if x not in DNA_bases]
# assign start values
input_num = 0
base_change = 0
base_change_seq_number = 0
base_change_seq_type = 0
base_change_all = 0
invalid_list = []
###################################
# parse input file
query = SeqIO.index(input_file, input_format)
for key, values in query.items():
    DNA_id = values.id
    DNA_description = values.description.replace(values.id, "").strip(" ?")
    DNA_sequence = str(values.seq.upper())
    #
    if cleaner_args == "Y":
        DNA_description_2 = str(DNA_description.replace(" ", "_"))
        for BadChr in string.punctuation:
            DNA_description_2 = DNA_description_2.replace(BadChr, "_")
            DNA_description_2 = re.compile(r'_{2,}').sub('_', DNA_description_2)
    #
    for base in Remove_Characters:
        if base in DNA_sequence:
            invalid_list.append(base)
            base_change_seq_number = str(DNA_sequence).count(base)
            DNA_sequence = DNA_sequence.replace(base, "N")
            base_change_seq_type += int(base_change_seq_number)
    #
    if verbosity_args == "Y":
        if cleaner_args == "Y":
            if output_format_args == 3:
                print(DNA_id + "_" + DNA_description_2 + "\t" + DNA_sequence)
            elif output_format_args == 1 or 2:
                print(">" + DNA_id + "_" + DNA_description_2 + "\n" + DNA_sequence)
        #        
        elif cleaner_args == "N":
            if output_format_args == 3:
                print(DNA_id + " " + DNA_description + "\t" + DNA_sequence)
            elif output_format_args == 1 or 2:
                print(">" + DNA_id + " " + DNA_description + "\n" + DNA_sequence)
    # cleaned sequences
    if cleaner_args == "Y":
        DNA_header = DNA_id + "_" + DNA_description_2
        final_seq = (SeqRecord(Seq(str(DNA_sequence)), id=DNA_header, description="", name=""))
    else:
        if output_format_args == 3:
            DNA_header = DNA_id + " " + DNA_description # 
            final_seq = (SeqRecord(Seq(str(DNA_sequence)), id=DNA_header, description="", name=""))
        elif output_format_args == 1 or 2:
            final_seq = (SeqRecord(Seq(str(DNA_sequence)), id=DNA_id, description=DNA_description, name=""))
    # write output
    SeqIO.write(final_seq, output_obj, output_format)
    #
    input_num += 1
#    
base_change_all += int(base_change_seq_type)
#
output_obj.close()
###################################
#
if not invalid_list:
    invalid_list = "none"
else:
    invalid_list = str(",".join(list(set(invalid_list))))
###################################
# Print info
print("#" * 60)
print("[Program]\t: " + program)
print("[Date]\t\t: "+ date_var_1.strftime("%Y-%m-%d %H:%M:%S"))
print("[Valid bases]\t: " + str(",".join(DNA_bases)))
print("[Input file]\t: "+ input_file)
print("\t\t  |_[Invalid bases]: " + invalid_list)
print("[Output file]\t: "+ output_file)
print("\t\t  |_[Output Seqs]: " + str(base_change_all) + " base change(s) in " +  str(input_num) + " sequence(s)")
print("[Header cleaned]: "+ clean_word)
print("#" * 60)
##################################
# end of script
