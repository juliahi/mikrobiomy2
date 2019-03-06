
from common import *
from common_parameters import *
from matplotlib import pyplot as plt
import numpy
import seaborn as sns
import pandas
import graph_stats
from parse_blast_output import count_qcov, count_qcov_len
import pysam

from plots import *
#from compare_heuristics import *


def compare_counts(graph, paths):
    counts = 0
    fcs = []
    for path in paths:
        c = get_path_count(graph, path)
        counts += c[0] + c[1]
        fcs.append(foldchange(*c))
    print "Sum = ", counts
    return fcs


def load_sga_scaffold(filename, sg, min_len=0):  # for my modified SGA scaffold results
    paths = []
    seqs = []
    not_in_graph = 0
    not_in_graph_longer = 0
    for name, seq in list(iter_fasta(filename)):
        path = name.split()[1][8:].split(';')[:-1]
        if sg.new_names is not None:
            path = [sg.new_names[x.strip('+-')] for x in path]
        else:
            path = [x.strip('+-') for x in path]

        if len(seq) > min_len:
            if len([1 for node in path if not sg.graph.has_node(node)]) > 0:
                not_in_graph += 1
                if len(path) > 1:
                    not_in_graph_longer += 1
            else:
                paths.append(path)
                seqs.append(seq)
    print "not loading %d paths that are not in graph any more. %d is not isolated" % \
          (not_in_graph, not_in_graph_longer)
    return paths, seqs


def print_stats(paths_list, seqs_list, counts_list, sums_list = None):
    print "No. paths: \tLength (bp): \tCounts0: \tCounts1: \tCounts: \t counts [%]: "
    norm = NPAIRS*2
    if sums_list is None:
        for paths, seqs, counts in zip(paths_list, seqs_list, counts_list):
            frac = 1. * sum(map(lambda x: sum(x), counts)) / norm * 100
            print '{:0,d}\t\t{:0,d}\t{:0,.1f}\t{:0,.1f}\t{:0,.1f}\t{:0,.2f}'.format(len(seqs), sum(map(len, seqs)),
                                                                                    sum(map(lambda x: x[0], counts)),
                                                                                    sum(map(lambda x: x[1], counts)),
                                                                                    sum(map(lambda x: sum(x), counts)),
                                                                                    frac)

    else:
        for paths, seqs, counts, sums in zip(paths_list, seqs_list, counts_list, sums_list):
            counts = [[i[0]*sums[0] / 1000000, i[1]*sums[1] / 1000000] for i in counts]
            frac = 1. * sum(map(lambda x: sum(x), counts)) / norm * 100
            print '{:0,d}\t\t{:0,d}\t{:0,.1f}\t{:0,.1f}\t{:0,.1f}\t{:0,.2f}'.format(len(seqs), sum(map(len, seqs)),
                                                                                    sum(map(lambda x: x[0], counts)),
                                                                                    sum(map(lambda x: x[1], counts)),
                                                                                    sum(map(lambda x: sum(x), counts)),
                                                                                    frac)


def write_unique(path_sets, seq_sets, names, prefix):
    counts = {}
    for seqs in seq_sets:
        for seq in seqs:
            if seq in counts: counts[seq] += 1
            else: counts[seq] = 1

    for paths, seqs, name in zip(path_sets, seq_sets, names):
        paths_u, seqs_u = [], []
        for path, seq in zip(paths, seqs):
            if counts[seq] == 1:
                paths_u.append(path)
                seqs_u.append(seq)
        print "Saving %d unique sequences in %s" % (len(paths_u),  prefix + '_' + name + '.fa')
        write_to_fasta(paths_u, seqs_u, prefix + '_' + name + '.fa')


##################################################

def check_read(r):
    if r.is_secondary: return False
    if r.is_unmapped: return False
    if r.is_supplementary: return False
    return True


def normalize(counts):
    s0 = 1.*sum([c[0] for c in counts])
    s1 = 1.*sum([c[1] for c in counts])
    print s0, s1
    result = [[c[0]*1000000/s0, c[1]*1000000/s1] for c in counts]
    return result, (s0, s1)


def get_sequences(filename, min_len):
    names = []
    seqs = []
    for name, seq in iter_fasta(filename):
        if len(seq) >= min_len:
            names.append(name)
            seqs.append(seq)
        else: print name, len(seq)
    return names, seqs


def filter_by_fc(names, seqs, counts, min_fc=MIN_FC):
    fnames = []
    fseqs = []
    fcounts = []
    for name, seq, count in zip(names, seqs, counts):
        if foldchange_compare(count[0], count[1], min_fc):
            fnames.append(name)
            fseqs.append(seq)
            fcounts.append(count)
    return fnames, fseqs, fcounts


def init_assembly(filename, bam_files):
    print filename
    names, seqs = get_sequences(filename, min_len=MIN_LENGTH)
    counts = load_coverage_from_bams(names, bam_files, conditions)
    print_stats([names], [seqs], [counts])
    # return names, seqs, counts, []
    counts, sums = normalize(counts)
    names_filter, seqs_filter, counts_filter = filter_by_fc(names, seqs, counts, min_fc=MIN_FC)
    print_stats([names_filter], [seqs_filter], [counts_filter], [sums])
    return names_filter, seqs_filter, counts_filter, sums


def load_coverage_from_bams(names, bam_files, conditions):
    covs = [[0, 0] for i in names]
    n = 0
    n2 = 0
    suppls = 0
    for filename in bam_files:
        with pysam.AlignmentFile(filename, 'rb') as bamfile:
            filereads = 0
            #print list(bamfile.get_index_statistics())[-5:]
            #print filename
            n += sum([x[1] for x in bamfile.get_index_statistics()]) #+ sum([x[2] for x in bamfile.get_index_statistics()])
            #chroms = set([x[0] for x in bamfile.get_index_statistics()])
            #print sum([x[1] for x in bamfile.get_index_statistics()]), sum([x[2] for x in bamfile.get_index_statistics()])
            #print "nocoordinate", bamfile.nocoordinate
            print "mapped & unmapped", bamfile.mapped, bamfile.unmapped, bamfile.mapped+ bamfile.unmapped
            print "chromosomes", len([x[1] for x in bamfile.get_index_statistics()])

            #names = [n.split()[0] for n in names]

            cond = conditions[filename.split('/')[-1][:15]]
            for i, name in enumerate(names):
                c = bamfile.count(name.split()[0], read_callback=check_read)
                #suppls += len([r for r in bamfile.fetch() if r.is_supplementary])
                covs[i][cond] += c
                filereads += c
                #if name.split()[0] in chroms:
                #    chroms.remove(name.split()[0])
            n2 += filereads
            #print filereads, float(filereads)/(bamfile.mapped+ bamfile.unmapped), float(filereads)/(NPAIRS*2)
            #print [covs[i][cond] for i in range(5)]

    print "coverage by bamfile.mapped", float(n)/(NPAIRS*2)
    print "coverage by iterating names", float(n2)/(NPAIRS*2)

    return covs


######## plot some properties #############

def plot_identical_seqs(seq_list, labels):
    seq_sets = [set(seqs) for seqs in seq_list]
    for idx1, idx2 in index_pairs(len(seq_sets)):
        ident = len(seq_sets[idx1].intersection(seq_sets[idx2]))
        print labels[idx1], labels[idx2], "Identical:", ident, \
            "S1 percent:", ident * 1.0 / len(seq_sets[idx1]), "S2 percent:", ident * 1.0 / len(seq_sets[idx2])

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    plot_venn(seq_sets, labels, ax, colors, 'Returned paths')


######## ploting blast results ############

def plot_blast_results(names, prefix, colors):
    N=1000
    qcovs, nhits = count_qcov([OUTPUTDIR + prefix + '/blast_' + name + '_random%d.tsv'%N for name in names], N)
    nhits = [[q if q < 10 else 10 for q in nhits1] for nhits1 in nhits]

    fig, ax = plt.subplots(1, 2, figsize=(16, 6))

    ax[0].hist(qcovs, label=names, color=[colors[x] for x in names])
    ax[0].legend()
    ax[0].set_xlabel("query coverage [%]")
    # ax[0,1].hist(qcovs, cumulative=True, normed=True, histtype='step', bins=50, label=names, linewidth=2)
    # ax[0,1].legend()

    ax[1].hist(nhits, label=names, color=[colors[x] for x in names])
    ax[1].set_xlabel("number of BLAST hits")
    ax[1].legend()


def all_blast_longest(names, prefix, colors):
        qcovs, nhits = count_qcov(['blast_' + prefix + '_' + name + 'longest1000.tsv' for name in names], 1000)
        nhits = [[q if q < 10 else 10 for q in nhits1] for nhits1 in nhits]

        fig, ax = plt.subplots(1, 2, figsize=(16, 6))

        ax[0].hist(qcovs, label=names, color=colors)
        ax[0].legend()
        ax[0].set_xlabel("query coverage [%]")
        # ax[0,1].hist(qcovs, cumulative=True, normed=True, histtype='step', bins=50, label=names, linewidth=2)
        # ax[0,1].legend()

        ax[1].hist(nhits, label=names, color=colors)
        ax[1].set_xlabel("number of BLAST hits")
        ax[1].legend()


def all_blast_length(labels, prefix, colors):
    df = count_qcov_len(['blast_' + prefix + '_' + name + '.tsv' for name in labels], labels)
    fig, ax = plt.subplots(1, 1, figsize=(16, 8))
    print labels

    sns.scatterplot(x="length", y="qcov", hue="label",
                    data=df, ax=ax, color=colors, #label=labels,
                    alpha=.7
                    )
    #sns.distplot(df, x="qcov", label=labels, color=colors, ax[1])
    #ax[1].legend()
    #ax[1].set_xlabel("query coverage [%]")
    #ax[2].hist(df["length"], label=labels, color=colors)
    #ax[2].legend()
    #ax[2].set_xlabel("query length [%]")


