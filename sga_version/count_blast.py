import pysam
import argparse
import matplotlib.pyplot as plt
import re
import pandas


def count_qcov(blast_filenames, N):
    qcovs = []
    nhits = []

    if type(N) == int:
        Ns = [N]*len(blast_filenames)
    else: Ns = N

    for filename, N in zip(blast_filenames, Ns):
        with open(filename, 'r') as blastfile:
            reads = {}
            reads_lengths = {}
            for line in blastfile:
                values = re.split(r'\s+', line)
            
                readname = values[0]
                coverage = float(values[12])

                if readname in reads:
                    reads[readname] += 1
                    reads_lengths[readname] = max(reads_lengths[readname], coverage)
                else:
                    reads[readname] = 1
                    reads_lengths[readname] = coverage
            zeros = [0 for x in xrange(N-len(reads_lengths))]
            qcovs.append(reads_lengths.values() + zeros)
            nhits.append(reads.values() + zeros)
    return qcovs, nhits


def count_qcov_len(blast_filenames, labels):
    df = pandas.DataFrame({'label': [], 'qcov': [], "length": []})

    for filename, label in zip(blast_filenames, labels):
        with open(filename, 'r') as blastfile:
            reads = {}
            reads_covs = {}
            for line in blastfile:
                values = re.split(r'\s+', line)

                readname = values[0]
                length = float(values[3])
                coverage = float(values[12])

                if readname in reads:
                    reads_covs[readname] = max(reads_covs[readname], coverage)
                else:
                    reads[readname] = length
                    reads_covs[readname] = coverage
            names = reads.keys()

            df = df.append(pandas.DataFrame({'label': [label]*len(names),
                                             'qcov': [reads_covs[r] for r in names],
                                             "length": [reads[r] for r in names]}
                                            ), ignore_index=True)
    return df


def main():
    parser = argparse.ArgumentParser(description='Count reads mapped with blast')
    parser.add_argument('regions', type=str, nargs='+',
                   help='blast files')
    parser.add_argument('-o', '--out', type=str, 
                   help='output file')
    args = parser.parse_args()

    count(args.regions, args.out)


if __name__ == "__main__":
   main()
