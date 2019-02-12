
from common import *
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import numpy
import pandas
import common
from load_filter_graph import give_time
import glob
import seaborn as sns

import pysam



def stats(seqs, counts):
    print "No. paths: \tLength (bp): \tCounts0: \tCounts1: \tCounts: \t % counts: "
    print len(seqs), "\t\t", \
        "\t", sum(map(len, seqs)), \
        "\t", sum(map(lambda x: x[0], counts)), \
        "\t", sum(map(lambda x: x[1], counts)), \
        "\t", sum(map(lambda x: sum(x), counts)), \
        "\t", 1.*sum(map(lambda x: sum(x), counts))/(2*NPAIRS)


def stats_after_norm(seqs, counts, sums):
    NPAIRS=90582474

    #reverse normalizing
    r0 = sum(map(lambda x: x[0], counts))*sums[0]/1000000
    r1 = sum(map(lambda x: x[1], counts))*sums[1]/1000000
    print "No. paths: \tLength (bp): \tCounts0: \tCounts1: \tCounts: \t % counts: \t normed counts[Mbp]"
    print len(seqs), "\t\t", \
        "\t", sum(map(len, seqs)), \
        "\t", r0, "\t", r1, "\t", r0+r1, \
        "\t", (r0+r1) / (2*NPAIRS), "\t", sum(map(lambda x: sum(x), counts))


def get_sequences(filename, min_len):
    names = []
    seqs = []
    for name, seq in iter_fasta(filename):
        if len(seq) >= min_len:
            names.append(name)
            seqs.append(seq)
        else: print name, len(seq)
    return names, seqs


def init_assembly(filename, bam_files, conditions, min_len, min_fc):
    print filename
    names, seqs = get_sequences(filename, min_len=min_len)
    counts = load_coverage_from_bams(names, bam_files, conditions)
    stats(seqs, counts)
    # return names, seqs, counts, []
    counts, sums = normalize(counts)
    names_filter, seqs_filter, counts_filter = filter_by_fc(names, seqs, counts, min_fc=min_fc)
    stats_after_norm(seqs_filter, counts_filter, sums)
    return names_filter, seqs_filter, counts_filter, sums


def check_read(r):
    if r.is_secondary: return False
    if r.is_unmapped: return False
    if r.is_supplementary: return False
    return True


def load_coverage_from_bams(names, bam_files, conditions):
    covs = [[0, 0] for i in names]
    n = 0
    n2 = 0
    suppls = 0
    for filename in bam_files:
        with pysam.AlignmentFile(filename, 'rb') as bamfile:
            filereads = 0
            print list(bamfile.get_index_statistics())[-5:]
            print filename
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
            print filereads, float(filereads)/(bamfile.mapped+ bamfile.unmapped), float(filereads)/(NPAIRS*2)
            #print [covs[i][cond] for i in range(5)]

    print "coverage by bamfile.mapped", float(n)/(NPAIRS*2)
    print "coverage by iterating names", float(n2)/(NPAIRS*2)

    return covs


def normalize(counts):
    s0 = 1.*sum([c[0] for c in counts])
    s1 = 1.*sum([c[1] for c in counts])
    result = [[c[0]*1000000/s0, c[1]*1000000/s1] for c in counts]
    return result, (s0, s1)


def load_base_coverage_from_bams(names, bam_files, conditions):
    covs = [[[], []] for _ in names]
    print len(names)
    for filename in bam_files:
        with pysam.AlignmentFile(filename, 'rb') as bamfile:
            print filename
            cond = conditions[filename.split('/')[-1][:15]]
            for i, name in enumerate(names):
                if covs[i][cond] is []:
                    covs[i][cond] = [c.nsegments
                                     for c in bamfile.pileup(name.split()[0], ignore_orphans=False)]
                else:
                    covs[i][cond] = [x + c.nsegments
                                     for x, c in zip(covs[i][cond],
                                                     bamfile.pileup(name.split()[0], ignore_orphans=False))]
    return covs


def coverage_variances(covs_lists, labels):
    heuristic = []
    variance0 = []
    dispersion0 = []
    variance1 = []
    dispersion1 = []

    for covs_list, name in zip(covs_lists, labels):
        print "counting ", name
        for cov0, cov1 in covs_list:
            heuristic += [name]
            variance0.append(variance(cov0))
            variance1.append(variance(cov1))
            dispersion0.append(dispersion(cov0))
            dispersion1.append(dispersion(cov1))

    df = pandas.DataFrame(data={'heuristic': heuristic, "variance0": variance0, "variance1": variance1,
                                "dispersion0": dispersion0, "dispersion1": dispersion1})
    return df


def coverage_variances_combine(names, bam_files, conditions, label):
    def coverage(seq_name):
        cov = []
        i = 0
        for c in bamfile.pileup(seq_name.split()[0], truncate=True, ignore_orphans=False):
            while i < c.reference_pos:
                cov.append(0)
                i += 1
            cov.append(c.nsegments)
        return cov

    covs = [[[], []] for _ in names]

    heuristic = []
    variance0 = []
    dispersion0 = []
    variance1 = []
    dispersion1 = []

    print len(names)
    for filename in bam_files:
        with pysam.AlignmentFile(filename, 'rb') as bamfile:
            print filename
            cond = conditions[filename.split('/')[-1][:15]]
            for i, name in enumerate(names[:5]):
                if len(covs[i][cond]) == 0:
                    covs[i][cond] = coverage(name)
                else:
                    new_c = coverage(name)
                    #if len(covs[i][cond]) != len(new_c):
                    #    print len(covs[i][cond]), len(new_c)
                    while len(covs[i][cond]) < len(new_c): covs[i][cond].append(0)
                    while len(covs[i][cond]) > len(new_c): new_c.append(0)
                    covs[i][cond] = [x + c for x, c in zip(covs[i][cond], new_c)]

    print len(covs)

    for cov0, cov1 in covs:
            heuristic += [label]
            variance0.append(variance(cov0))
            variance1.append(variance(cov1))
            dispersion0.append(dispersion(cov0))
            dispersion1.append(dispersion(cov1))

    df = pandas.DataFrame(data={'heuristic': heuristic, "variance0": variance0, "variance1": variance1,
                                "dispersion0": dispersion0, "dispersion1": dispersion1})
    print df.head()
    return df


def filter_by_fc(names, seqs, counts, min_fc):
    fnames = []
    fseqs = []
    fcounts = []
    for name, seq, count in zip(names, seqs, counts):
        if foldchange_compare(count[0], count[1], min_fc):
            fnames.append(name)
            fseqs.append(seq)
            fcounts.append(count)
    return fnames, fseqs, fcounts


def init_megahit(conds, min_len, min_fc):
    # dir = "/mnt/chr4/mikrobiomy-2/megahit_results/all/"
    dir = "/home/julia/mikrobiomy_results/megahit/"
    filename = dir + "long_contigs_200.fa"
    bam_files = glob.glob(dir + "*sorted.bam")
    return init_assembly(filename, bam_files, conds, min_len, min_fc)


def foldchange_distr(count_lists, labels, colors):
    fig, ax = plt.subplots(1, figsize=(15, 9))
    for i in xrange(len(count_lists)):
        values = map(lambda count: math.log(foldchange(count[0], count[1]), 2), count_lists[i])
        #values = map(lambda count: foldchange(count[0], count[1]), count_lists[i])
        print "%s (len %d): '-Inf' values: %f, +Inf: %f" % (labels[i], len(count_lists[i]),
                                                            1.*len([x for x in values if x <= -12])/len(values),
                                                            1.*len([x for x in values if x >= 12])/len(values))

        sns.distplot([x for x in values if -12 < x < 12], rug=False, hist=False, kde=True,
                     label=labels[i], ax=ax, color=colors[i])

    ax.set(xlabel="log2(fold change)")
    ax.set_xlim(-12, 12)
    plt.legend()
    return ax



