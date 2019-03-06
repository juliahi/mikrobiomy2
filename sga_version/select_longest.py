# iter over fasta sequences

from itertools import groupby
import random
import sys


def fasta_iter(fasta_file):
    for header, group in groupby(fasta_file, lambda line: line[0] == ">"):
        if header:
            line = group.next()
            ensembl_id = line[1:].strip()
        else:
            sequence = ''.join(line.strip() for line in group)
            yield ensembl_id, sequence


def filter_file(inputf, outputf, N):
    sequences = []
    with open(outputf, 'w+') as f:
        for name, seq in fasta_iter(open(inputf)):
            sequences.append((name, seq))

        for name, seq in sorted(sequences, key=lambda x: len(x[1]), reverse=True)[:N]:
            f.write(">%s\n%s\n" % (name, seq))


filter_file(sys.argv[1], sys.argv[2], int(sys.argv[3]))


