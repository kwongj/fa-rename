#!/usr/bin/env python3

# Usage
import argparse
from argparse import RawTextHelpFormatter
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO
import sys

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
	usage='\n  %(prog)s [--ids new_names.txt] FASTA > new.fasta')
parser.add_argument('fasta', metavar='FASTA', nargs=1, help='original FASTA file')
parser.add_argument('--ids', metavar='FILE', required=True, nargs=1, help='specify two column tab-separated file with [oldnames] [newnames]')
parser.add_argument('--out', metavar='FILE', help='specify output file (default = stdout)')
parser.add_argument('--version', action='version', version='%(prog)s v0.1')
args = parser.parse_args()

# Rename leaf nodes
newseqs = []
new_names = tab2dict(args.ids[0], '\t')
for record in SeqIO.parse(args.fasta[0], 'fasta'):
	newid = new_names[record.id]
	newseqs.append(SeqRecord(record.seq, id=newid, description=''))

# Write masked alignment to file or print to stdout
if args.out:
	msg('Masked sequences saved to "{}" ... '.format(outfile))
	SeqIO.write(newseqs, outfile, 'fasta')
else:
	seqFILE = StringIO.StringIO()
	SeqIO.write(newseqs, seqFILE, 'fasta')
	output = seqFILE.getvalue().rstrip()
	print(output)

sys.exit(0)
