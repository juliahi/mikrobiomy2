
from common import *
from common_parameters import *
from SGAFilter import *
import numpy
import networkx as nx
import seaborn as sns
import pandas

from graph_stats import get_path_count, get_path_coverage
from plots import *

import cPickle


def stats(graph, paths, seqs):
    counts = map(lambda x: get_path_count(graph, x), paths)

    unnormed = (sum(map(lambda x: x[0], counts)) * graph.sums[0] + \
           sum(map(lambda x: x[1], counts)) * graph.sums[1]) / 1000000
    print "No. paths: \tNo. vertices: \tLength (bp): \tCounts0: \tCounts1: \tCounts:"
    print len(paths), "\t\t", sum([len(path) for path in paths]), \
        "\t", sum(map(len, seqs)), \
        "\t", sum(map(lambda x: x[0], counts)), \
        "\t", sum(map(lambda x: x[1], counts)), \
        "\t", sum(map(lambda x: sum(x), counts))
    print "\t Counts not normed \t Counts not normed [fraction]\t counts from graph [fraction]"
    print unnormed, "\t", unnormed / (2*NPAIRS), \
        "\t", sum(map(lambda x: sum(x), counts)) / sum(graph.node_nreads_sums())


def save_paths_to_fasta(sg, paths, seqs, out_file):        # writing to FASTA
        counts = 0.
        length = 0
        total = 0
        with open(out_file, 'w+') as f:
            for path, seq in zip(paths, seqs):
                counts_tmp = get_path_count(sg, path)
                counts += counts_tmp[0] + counts_tmp[1]
                length += len(seq)
                total += 1
                f.write('>%s counts1=%f counts2=%f foldchange=%f\n%s\n' % (
                        PATH_SEP.join(path), counts_tmp[0], counts_tmp[1],
                        foldchange(counts_tmp[0], counts_tmp[1]), seq))
        print "Saved %d sequences with %d counts of total length %d" % (total, counts, length)


def load_info_from_fasta(filename):
    paths = []
    seqs = []
    counts = []

    for name, seq in iter_fasta(filename):
        name, c1, c2, fc = name.split()
        paths.append(name.split(PATH_SEP))
        counts.append((float(c1.split('=')[1]), float(c2.split('=')[1])))
        seqs.append(seq)
    return paths, seqs, counts


def load_paths_from_fasta(filename):
    paths = []
    seqs = []

    for name, seq in iter_fasta(filename):
        name = name.strip().split()
        paths.append(name.split(PATH_SEP))
        seqs.append(seq)
    return paths, seqs


### Filtering selected paths ###

def filter_by_fc(graph, paths, sequences, min_fc, new_names=False):
        if len(paths) == 0:
            return []

        counts = 0
        filtered_paths = []
        filtered_seqs = []
        for path, s in zip(paths, sequences):
            counts_tmp = get_path_count(graph, path, new_names)
            if not foldchange_compare(counts_tmp[0], counts_tmp[1], min_fc):
                continue
            counts += counts_tmp[0] + counts_tmp[1]
            filtered_paths.append(path)
            filtered_seqs.append(s)

        print "Leave %d sequences with %d counts" % (len(filtered_paths), counts)
        return filtered_paths, filtered_seqs


def filter_by_seqlen(paths, seqs, min_len):
    filtered_paths = []
    filtered_seqs = []
    for p, s in zip(paths, seqs):
        if len(s) >= min_len:
            filtered_paths.append(p)
            filtered_seqs.append(s)
    return filtered_paths, filtered_seqs


def filter_by_pathlength(paths, seqs, l=1):
    filtered_paths = []
    filtered_seqs = []
    for p, s in zip(paths, seqs):
        if len(p) > l:
            filtered_paths.append(p)
            filtered_seqs.append(s)
    return filtered_paths, filtered_seqs


def remove_identical(path_sets, seq_sets, count_sets):      # remove common 1-node paths
    sets = []
    for paths in path_sets:
        s = set([x[0] for x in paths if len(x) == 1])
        sets.append(s)
        print "1-node paths:", len(s)
    int = set.intersection(*sets)
    print "Intersection:", len(int)

    paths2 = []
    seqs2 = []
    counts2 = []
    for paths, seqs, counts in zip(path_sets, seq_sets, count_sets):
        p2 = []
        s2 = []
        c2 = []
        for p, s, c in zip(paths, seqs, counts):
            if len(p) > 1 or p[0] not in int:
                p2.append(p)
                s2.append(s)
                c2.append(c)
        paths2.append(p2)
        seqs2.append(s2)
        counts2.append(c2)

    return paths2, seqs2, counts2


### Comparing used vertices and edges, pairs of edges ###

def plot_used_vertices(paths, names, colors):
    node_sets = [set([node for path in paths_set for node in path]) for paths_set in paths]
    for idx1, idx2 in index_pairs(len(paths)):
        print names[idx1], names[idx2]
        jaccard_index(node_sets[idx1], node_sets[idx2])

    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    plot_venn(node_sets[:3], names[:3], ax[0], colors, 'Comparison of used vertices\n')
    plot_venn(node_sets[1:], names[1:], ax[1], colors, 'Comparison of used vertices\n')


def plot_used_edges(path_sets, names, colors):
    def get_edges_from_list(paths):
        all_edges = []
        for names in paths:
            edges = cons_pairs(names)
            all_edges += edges
        return set(all_edges)

    edge_sets = [set(get_edges_from_list(paths)) for paths in path_sets]
    for idx1, idx2 in index_pairs(len(edge_sets)):
        print names[idx1], names[idx2]
        jaccard_index(edge_sets[idx1], edge_sets[idx2])
    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    plot_venn(edge_sets[:3], names[:3], ax[0], colors, 'Comparison of used edges\n')
    plot_venn(edge_sets[1:], names[1:], ax[1], colors, 'Comparison of used edges\n')


def plot_used_pairs(path_sets, names, colors):
    def get_pairs_of_edges(paths):
        all_pairs = []
        for path in paths:
            if len(path) >= 3:
                all_pairs += [(path[i], path[i + 1], path[i + 2]) for i in xrange(len(path) - 2)]
        return set(all_pairs)

    edge_sets = [set(get_pairs_of_edges(paths)) for paths in path_sets]
    for idx1, idx2 in index_pairs(len(edge_sets)):
        print names[idx1], names[idx2]
        jaccard_index(edge_sets[idx1], edge_sets[idx2])
    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    plot_venn(edge_sets[:3], names[:3], ax[0], colors, 'Comparison of used pairs of edges\n')
    plot_venn(edge_sets[1:], names[1:], ax[1], colors, 'Comparison of used pairs of edges\n')


def plot_identical_paths(path_sets, names, colors):
    path_sets = [set(map(tuple, paths)) for paths in path_sets]
    for idx1, idx2 in index_pairs(len(path_sets)):
        ident = len(path_sets[idx1].intersection(path_sets[idx2]))
        print names[idx1], names[idx2], "Identical:", ident, \
            "S1 percent:", ident * 1.0 / len(path_sets[idx1]), "S2 percent:", ident * 1.0 / len(path_sets[idx2])

    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    plot_venn(path_sets[:3], names[:3], ax[0], colors, 'Returned paths')
    plot_venn(path_sets[1:4], names[1:4], ax[1], colors, 'Returned paths')


### Plot fold-changes ###

def plot_foldchanges(count_lists, labels, colors):
    fig, ax = plt.subplots(1, figsize=(15, 9))
    for i in xrange(len(count_lists)):
        values = count_lists[i]
        print "%s: '-Inf': %.4f, +Inf: %.4f" % (labels[i], 1.*len([x for x in values if x[0] == 0])/len(values),
                                                       1.*len([x for x in values if x[1] == 0])/len(values))
        values = [math.log(foldchange(*c), 2) for c in count_lists[i] if c[0] != 0 and c[1] != 0]
        sns.distplot(values, rug=False, hist=False, kde=True,
                     label=labels[i], ax=ax, color=colors[labels[i]])

    ax.set(xlabel="log2(fold change)")
    #ax.set_xlim(-15,15)
    plt.legend()
    # df = pandas.DataFrame(data={'heuristic': h, 'lfc': lfcs})
    # return df
    return ax


### Plot node coverage ###

def prepare_path_coverage_data(sg, path_sets, names, min_length=1):
    heuristic = []
    variance0 = []
    dispersion0 = []
    variance1 = []
    dispersion1 = []
    c=0

    for paths, name in zip(path_sets, names):
        print "counting ", name
        for path in paths:
            if len(path) >= min_length:
                heuristic += [name]
                cov0, cov1 = get_path_coverage(sg, path)
                variance0.append(variance(cov0))
                variance1.append(variance(cov1))
                dispersion0.append(dispersion(cov0))
                dispersion1.append(dispersion(cov1))
                c+=1
                if c % 10000 == 0: print c

    df = pandas.DataFrame(data={'heuristic': heuristic, "variance0": variance0, "variance1": variance1,
                                "dispersion0": dispersion0, "dispersion1": dispersion1})
    return df


def plot_path_coverage_variance(df, colors, remove=0.05):
    l = len(df['heuristic'])
    df2 = pandas.DataFrame(data={'heuristic': list(df['heuristic']) + list(df['heuristic']),
                                 "variance": list(df['variance0']) + list(df['variance1']),
                                 "dispersion": list(df['dispersion0']) + list(df["dispersion1"]),
                                 "condition": [0]*l + [1]*l})

    df3 = df2.copy()
    for name in df3.heuristic.unique():
        qv0 = df3[(df3.condition == 0) & (df3.heuristic == name)]["variance"].quantile(1-remove)
        qv1 = df3[(df3.condition == 1) & (df3.heuristic == name)]["variance"].quantile(1-remove)
        print name, qv0, qv1
        df3 = df3[(df3.heuristic != name) | (((df3.condition == 0) & (df3.variance < qv0))
                  | ((df3.condition == 1) & (df3.variance < qv1)))]

    fig, ax = plt.subplots(2, 1, figsize=(20, 18))
    sns.violinplot(x="heuristic", y="variance", hue="condition", data=df3,
                   ax=ax[0],  bw=.2, cut=0)

    df3 = df2
    df3[df3.dispersion.isna()]["dispersion"] = 0

    for name in df3.heuristic.unique():
        qd0 = df3[(df3.condition == 0) & (df3.heuristic == name)]["dispersion"].quantile(1-remove)
        qd1 = df3[(df3.condition == 1) & (df3.heuristic == name)]["dispersion"].quantile(1-remove)
        print name, qd0, qd1

        df3 = df3[(df3.heuristic != name) | (((df3.condition == 0) & (df3.dispersion < qd0))
                                             | ((df3.condition == 1) & (df3.dispersion < qd1)))]
    sns.violinplot(x="heuristic", y="dispersion", hue="condition", data=df3,
                   ax=ax[1], bw=.2, cut=0)
    # ax[1].set(yscale='log')
    # ax[1].set(xlabel="variance of log2(fold change)")




# TODO

"""
def compare_counts(graph, paths):
    counts = 0
    fcs = []
    for path in paths:
        c = get_path_count(graph, path)
        counts += c[0] + c[1]
        fcs.append(foldchange(*c))
    print "Sum = ", counts
    return fcs


def lfc_variance(fc, total_fc):
    return sum([(math.log(x, 2)-math.log(total_fc, 2))**2 for x in fc])/len(fc)
    #return numpy.var([math.log(x, 2) for x in fc])
    #return numpy.var(fc)


def fc_weighted_variance(fc, weights):
    s = 1.*sum(weights)
    m = numpy.mean([math.log(x, 2)*y/s for x, y in zip(fc, weights)])
    return sum([y/s*((math.log(x, 2)*y/s - m)**2) for x, y in zip(fc, weights)])
    # return numpy.var([1.*x*y/s for x, y in zip(fc, weights)])


# Wrapping functions

def filter_islands(sg, min_len, min_fc):
    node_list = []
    sum_len = 0
    for node, deg in sg.graph.degree():
        if deg == 0:
            node_list.append(node)
            sum_len += sg.graph.node[node]["length"]
    fc_count = len([1 for node in node_list
                    if foldchange_compare(sg.counts[node][0], sg.counts[node][1], min_fc) and
                    sg.graph.node[node]["length"] > min_len])

    print "Remove %d islands of length %d, %d with |foldchange| >= %f and length >= %d" % \
          (len(node_list), sum_len, fc_count, min_fc, min_len)

    # removing nodes
    sg.remove_nodes(node_list)
    return sg


def get_nodes_foldchange(sg, paths, single=True):
    # single - include single-node paths
    if single:
        fcs = [map(sg.node_foldchange, nodes) for nodes in paths]
        lens = [map(lambda x:sg.graph.nodes[x]["length"], nodes) for nodes in paths]
    else:
        fcs = [map(sg.node_foldchange, nodes) for nodes in paths if len(nodes) > 1]
        lens = [map(lambda x:sg.graph.nodes[x]["length"], nodes) for nodes in paths if len(nodes) > 1]
    return fcs, lens

def prepare_variance_data(sg, path_sets, names):
    h = []
    variance = []
    weighted = []
    dispersion = []
    dispersion_fc = []

    for paths, name in zip(path_sets, names):
        paths = [p for p in paths if len(p) > 0]
        h += [name] * len(paths)
        fcs, lengths = get_nodes_foldchange(sg, paths)
        path_fc = map(lambda path: foldchange(*get_path_count(sg, path)), paths)
        lfc_var = [lfc_variance(*x) for x in zip(fcs, path_fc)]
        wfc = map(lambda x: fc_weighted_variance(x[0], x[1]), zip(fcs, lengths))
        variance += lfc_var
        weighted += wfc
        dispersion += map(lambda x: x/numpy.mean(lfc_var), lfc_var)
        dispersion_fc += [x/math.log(y, 2) if y != 1 else None for x, y in zip(lfc_var, path_fc)]
    df = pandas.DataFrame(data={'heuristic': h, "variance": variance, "weighted": weighted,
                                "dispersion_mean": dispersion, "dispersion_real": dispersion_fc})
    return df


def all_fc_variance(df, names, colors):
    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    sns.violinplot(x="heuristic", y="variance", hue="heuristic", data=df[df["variance"] > 0.01],
                            palette="muted", ax=ax[0], color=colors)
    # ax[0].set(yscale='log')
    # ax[0].set(xlabel="variance of log2(fold change)")

    sns.violinplot(x="heuristic", y="weighted", hue="heuristic", data=df[df["variance"] > 0.01],
                            palette="muted", ax=ax[1], color=colors)
    ax[1].set(yscale='log')
    # ax[1].set(xlabel="variance of log2(fold change)")



def coverage_heatmap(df, name):
    #TODO
    pass


def all_fc_dispersion(df, names, colors):
    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    sns.violinplot(x="heuristic", y="dispersion_mean", hue="heuristic", data=df[df["variance"] > 0.01],
                            palette="muted", ax=ax[0], color=colors)
    # ax[0].set(yscale='log')
    ax[0].set(xlabel="dispersion of log2(fold change)")

    sns.violinplot(x="heuristic", y="dispersion_real", hue="heuristic", data=df[df["variance"] > 0.01],
                            palette="muted", ax=ax[1], color=colors)
    # ax[1].set(yscale='log')



def check_pairing(path_sets, to_old_names):
    result = []
    for paths in path_sets:
        same_path = 0
        same_node = 0
        different_paths = 0
        reads = set()
        for path in paths:
            tmp = set()
            for node in path:
                tmp2 = set()
                for read in to_old_names[node].split(NODE_SEP):
                    if read[:-2] in tmp2:
                        same_node += 2
                        tmp2.remove(read[:-2])
                    elif read[:-2] in tmp:
                        same_path += 2
                        tmp.remove(read[:-2])
                    elif read[:-2] in reads:
                        different_paths += 2
                        reads.remove(read[:-2])
                    else:
                        tmp2.add(read[:-2])
                tmp |= tmp2
            reads |= tmp
        result.append((same_node, same_path, different_paths, len(reads)))
    return result
"""