import Bio
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Filter FASTA file using ids from a file')
parser.add_argument('-f', help='FASTA file', required=True)
parser.add_argument('-i', help='file with ids to keep', required=True)
parser.add_argument('-o', help='output FASTA file', required=True)

args = parser.parse_args()

## Read ids
ids = []
for line in open(args.i, 'r'):
    ids.append(line.rstrip())

## Scan FASTA file
outf = open(args.o, 'w')
for rec in SeqIO.parse(args.f,'fasta'):
    if rec.id in ids:
        SeqIO.write(rec, outf,"fasta")
outf.close()
