#!/usr/bin/env python3

# Usage
import argparse
from argparse import RawTextHelpFormatter
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO
import os
import sys
import re

# Functions
# Log a message to stderr
def msg(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
	msg(*args, **kwargs)
	sys.exit(1);

# Check file exists
def check_file(f):
	return os.path.isfile(f)

# Check if file is in FASTA format
def check_fasta(f):
	if not os.path.isfile(f) or os.path.getsize(f) < 1:
		return False
	with open(f, 'r') as fasta:
		if fasta.readline()[0] != '>':						# Check if header starts with ">"
			return False
		for line in fasta:
			line = line.strip()
			if not line or line[0] == '>':	
				continue
			if bool(re.search('[^ACTGactgNn?\-]', line)):	# Check if there are non-nucleotide characters in sequence
				return False
	return True

def tab2dict(tab, sep):
	dict = {}
	with open(tab, mode='r') as file:
		table = csv.reader(file, delimiter=sep)
		for row in table:
			dict[row[0]] = row[1]
	return dict

parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='Rename headers/sequence IDs in multi-FASTA file\n',
	usage='\n  %(prog)s [--tab new_names.txt] FASTA > new.fasta')
parser.add_argument('fasta', metavar='FASTA', nargs=1, help='original FASTA file')
parser.add_argument('--ids', metavar='FILE', required=True, nargs=1, help='specify tab-separated file with [oldnames] [newnames]')
parser.add_argument('--out', metavar='FILE', nargs=1, help='specify output file (default = stdout)')
parser.add_argument('--version', action='version', version='%(prog)s v0.1')
args = parser.parse_args()

# Check input/output files
if not check_file(args.fasta[0]):
	err('ERROR: Cannot find "{}". Check file exists in the specified directory.'.format(args.fasta[0]))
if not check_fasta(args.fasta[0]):
	err('ERROR: Check "{}" is in FASTA format.'.format(args.fasta[0]))
if not check_file(args.ids[0]):
	err('ERROR: Cannot find "{}". Check file exists in the specified directory.'.format(args.ids[0]))
if args.out:
	if check_file(args.out[0]):
		err('ERROR: "{}" already exists.'.format(args.out[0]))

# Rename leaf nodes
newseqs = []
new_names = tab2dict(args.ids[0], '\t')
for record in SeqIO.parse(args.fasta[0], 'fasta'):
	newid = new_names[record.id]
	newseqs.append(SeqRecord(record.seq, id=newid, description=''))

# Write masked alignment to file or print to stdout
if args.out:
	msg('Masked sequences saved to "{}" ... '.format(args.out[0]))
	SeqIO.write(newseqs, args.out[0], 'fasta')
else:
	seqFILE = StringIO.StringIO()
	SeqIO.write(newseqs, seqFILE, 'fasta')
	output = seqFILE.getvalue().rstrip()
	print(output)

sys.exit(0)

