
from common import *
from SGAFilter import *
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import numpy
import networkx as nx
import seaborn as sns
import pandas
import common
import graph_stats
from load_filter_graph import give_time
import heuristics
from count_blast import count_qcov, count_qcov_len

import cPickle


jaccard_header = "Set1 \ set2\tSet1 n Set2\tSet2 \ set1\tJaccard index"


def jaccard_index(set1, set2):
    wspolne = len(set1.intersection(set2))
    idx = 1. * wspolne / len(set1.union(set2))
    #print "Set1\set2=", len(set1) - wspolne, '\t', "Set1 n Set2=", wspolne, '\t', "Set2\set1=", len(set2) - wspolne
    #print "%n\t%d\t%d\t%f" % (len(set1) - wspolne, wspolne, len(set2) - wspolne, idx)
    print '{:0,d}\t{:0,d}\t{:0,d}\t{:0,.4f}'.format(len(set1) - wspolne, wspolne, len(set2) - wspolne, idx)
    return idx


def get_paths_from_file(filename, sg, min_len=0):
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
    print "not loading %d paths that are not in graph any more. %d is not isolated" % (not_in_graph, not_in_graph_longer)
    return paths, seqs


def get_edges_from_list(paths):
    all_edges = []
    for names in paths:
        edges = cons_pairs(names)
        all_edges += edges
    return set(all_edges)


# def compare_used_edges(file1, file2):
#     edges1 = get_edges_from_list(file1)
#     edges2 = get_edges_from_list(file2)
#
#     index = jaccard_index(edges1, edges2)
#     #print "#edges1=%d\t#edges2=%d\tJaccard index=%f" % (len(edges1), len(edges2), index)


def get_pairs_of_edges(paths):
    all_pairs = []
    for path in paths:
        if len(path) >= 3:
            all_pairs += [(path[i], path[i+1], path[i+2]) for i in xrange(len(path)-2)]
    return set(all_pairs)

#
# def compare_pairs_of_edges(paths1, paths2):
#     # edges1 = get_pairs_of_edges(paths1)
#     # edges2 = get_pairs_of_edges(paths2)
#     #
#     # index = jaccard_index(edges1, edges2)
#     # print "#pairs1=%d\t#pairs2=%d\tJaccard index=%f" % (len(edges1), len(edges2), index)


# def compare_identical(paths1, paths2):
#     s1 = set(map(tuple, paths1))
#     s2 = set(map(tuple, paths2))
#     ident = len(s1.intersection(s2))
#
#     print "Identical:", ident, "S1 percent:", ident*1.0/len(s1), "S2 percent:", ident*1.0/len(s2)


def compare_counts(graph, paths):
    counts = 0
    fcs = []
    for path in paths:
        c = get_path_count(graph, path)
        counts += c[0] + c[1]
        fcs.append(foldchange(*c))
    print "Sum = ", counts
    return fcs


def get_info_from_file(filename):
    info = []
    with open(filename) as f:
        while True:
            line = f.readline()
            if not line: break
            if line.startswith(">"):
                name, counts1, counts2, fc = line.split()[0]
                name = name[1:]
                length = length(f.readline().strip())
                info.append((name, counts1, counts2, fc, length))
    return info


def filter_lengths(paths, seqs, min_len):
    filtered_paths = []
    filtered_seqs = []
    for p, s in zip(paths, seqs):
        if len(s) >= min_len:
            filtered_paths.append(p)
            filtered_seqs.append(s)
    return filtered_paths, filtered_seqs


def filter_pathlength(paths, seqs, l=1):
    filtered_paths = []
    filtered_seqs = []
    for p, s in zip(paths, seqs):
        if len(p) > l:
            filtered_paths.append(p)
            filtered_seqs.append(s)
    return filtered_paths, filtered_seqs


def stats(graph, paths, seqs):
    # print "No. paths:", len(paths), "No. vertices:", sum([len(path) for path in paths]), \
    #     "Length (bp): ", sum(map(len, seqs)), \
    #     "Counts:", sum(map(lambda x:sum(get_path_count(graph, x)), paths))
    counts = map(lambda x: get_path_count(graph, x), paths)
    print "No. paths: \tNo. vertices: \tLength (bp): \tCounts0: \tCounts1: \tCounts:"
    print len(paths), "\t\t", sum([len(path) for path in paths]), \
        "\t", sum(map(len, seqs)), \
        "\t", sum(map(lambda x: x[0], counts)), \
        "\t", sum(map(lambda x: x[1], counts)), \
        "\t", sum(map(lambda x: sum(x), counts))


def get_nodes_foldchange(sg, paths, single=True):
    # single - include single-node paths
    if single:
        fcs = [map(sg.node_foldchange, nodes) for nodes in paths]
        lens = [map(lambda x:sg.graph.nodes[x]["length"], nodes) for nodes in paths]
    else:
        fcs = [map(sg.node_foldchange, nodes) for nodes in paths if len(nodes) > 1]
        lens = [map(lambda x:sg.graph.nodes[x]["length"], nodes) for nodes in paths if len(nodes) > 1]
    return fcs, lens


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

def init_and_load():
    #conds = {'6683_16-06-2015': 1, '6685_04-06-2015': 0, '6685_16-06-2015': 1,
    #         '6690_04-06-2015': 0, '6690_16-06-2015': 1, '6695_04-06-2015': 0,
    #         '6695_16-06-2015': 1, '6704_04-06-2015': 0, '6704_16-06-2015': 1}
    test_name = 'sga_test_full_notrim_paired_reversed'
    #folder = "/mnt/chr7/data/julia/" + test_name
    suf = ".preprocessed_qf5.ec.filter.pass.rmdup"
    #filename = folder + "/merged" + suf + "_31.asqg"
    sg = cPickle.load(open('simplified1_'+test_name+suf+'_sga200_rename.pickle', 'rb'))
    #sg.write_to_asqg('renamed1_200-graph.asqg', contigs='renamed1_200-contigs.fa', rename=True)

    graph_stats.short_summary(sg)
    n_components, _ = sg.connected_components()
    print "Number of components", n_components
    return sg


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


def init_longest(sg, min_len, min_fc, outname):
    # Heuristic 1: take longest
    longest = heuristics.take_longest(sg)
    longest_seqs = sg.save_path_sequences(longest, outname+'_longest.fa')
    longest_filter, longest_filter_seqs = filter_paths(sg, longest, longest_seqs, min_fc=min_fc)
    longest_filter200, longest_filter200_seqs = filter_lengths(longest_filter, longest_filter_seqs, min_len)

    print "max_len before FC filter\t %d" % len(longest)
    print "max_len after fc filter"
    stats(sg, longest_filter, longest_filter_seqs)
    print "max_len after lengths filter"
    stats(sg, longest_filter200, longest_filter200_seqs)

    return longest_filter200, longest_filter200_seqs


def init_longest_fc(sg, min_len, min_fc, outname):
    # Heuristic 2: take longest until foldchange
    longest_fc = heuristics.take_longest_minfc(sg, min_fc)
    longest_fc_seqs = sg.save_path_sequences(longest_fc, outname + '_longestfc.fa')
    longest_fc200, longest_fc200_seqs = filter_lengths(longest_fc, longest_fc_seqs, min_len)

    print "max_lenfc"
    stats(sg, longest_fc, longest_fc_seqs)
    print "max_lenfc after lengths filter"
    stats(sg, longest_fc200, longest_fc200_seqs)

    return longest_fc200, longest_fc200_seqs


def init_best_fc(sg, min_len, min_fc, outname):
    # Heuristic 3: take best foldchange until > FC
    best_fc = heuristics.take_best_fc(sg, min_fc)
    best_fc_seqs = sg.save_path_sequences(best_fc, outname + '_bestfc.fa')
    best_fc200, best_fc200_seqs = filter_lengths(best_fc, best_fc_seqs, min_len)

    print "best_fc"
    stats(sg, best_fc, best_fc_seqs)
    print "best_fc after lengths filter"
    stats(sg, best_fc200, best_fc200_seqs)

    return best_fc200, best_fc200_seqs


def init_sga(sg, filename, min_len, min_fc):
    sga_paths, sga_seqs = get_paths_from_file(filename, sg, min_len=min_len)

    print "sga"
    stats(sg, sga_paths, sga_seqs)
    sga_filter = filter_paths(sg, sga_paths, sga_seqs, min_fc=min_fc)
    print "sga after fc filter"
    stats(sg, sga_filter[0], sga_filter[1])
    return sga_filter


def mix_color(c1, c2, c3=None):
    if c3 is None:
        return (c1[0] + c2[0])/2, (c1[1] + c2[1])/2, (c1[2] + c2[2])/2
    return (c1[0] + c2[0] + c3[0]) / 3, (c1[1] + c2[1] + c3[1]) / 3, (c1[2] + c2[2] + c3[2]) / 3


def plot_venn(sets, names, ax, colors):
    v = venn3(sets, set_labels=names, ax=ax)
    #c = venn3_circles(sets, ax=ax)

    v.get_patch_by_id('100').set_color(colors[0])
    v.get_patch_by_id('010').set_color(colors[1])
    v.get_patch_by_id('001').set_color(colors[2])
    v.get_patch_by_id('100').set_alpha(0.9)
    v.get_patch_by_id('010').set_alpha(0.9)
    v.get_patch_by_id('001').set_alpha(0.9)
    try:
        v.get_patch_by_id('110').set_color(mix_color(colors[0], colors[1]))
        v.get_patch_by_id('011').set_color(mix_color(colors[1], colors[2]))
        v.get_patch_by_id('101').set_color(mix_color(colors[0], colors[2]))
        v.get_patch_by_id('110').set_alpha(1)
        v.get_patch_by_id('011').set_alpha(1)
        v.get_patch_by_id('101').set_alpha(1)
        v.get_patch_by_id('111').set_color(mix_color(colors[0], colors[1], colors[2]))
        v.get_patch_by_id('111').set_alpha(0.9)
    except Exception:
        pass


def vertices_venn(sets, names, ax, colors):
    plot_venn(sets, names, ax, colors)
    ax.set_title('Comparison of used vertices\n')


def identical_venn(sets, names, ax, colors):
    plot_venn(sets, names, ax, colors)
    ax.set_title('Comparison of identical returned paths\n')


def all_vertices(paths, names, colors):
    node_sets = [set([node for path in paths_set for node in path]) for paths_set in paths]
    for idx1, idx2 in index_pairs(len(paths)):
        print names[idx1], names[idx2]
        #compare_used_vertices(paths[idx1], paths[idx2])
        jaccard_index(node_sets[idx1], node_sets[idx2])

    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    vertices_venn(node_sets[:3], names[:3], ax[0], colors[:3])
    vertices_venn(node_sets[1:], names[1:], ax[1], colors[1:])


def remove_identical(path_sets, seq_sets):
    sets = []
    for paths in path_sets:
        s = set([x[0] for x in paths if len(x) == 1])
        sets.append(s)
        print "1-node paths:", len(s)
    int = set.intersection(*sets)
    print "Intersection:", len(int)

    paths2 = []
    seqs2 = []
    for paths, seqs in zip(path_sets, seq_sets):
        p2 = []
        s2 = []
        for p, s in zip(paths, seqs):
            if len(p) > 1 or p[0] not in int:
                p2.append(p)
                s2.append(s)
        paths2.append(p2)
        seqs2.append(s2)

    return paths2, seqs2


def all_edges(path_sets, names, colors):
    edge_sets = [set(get_edges_from_list(paths)) for paths in path_sets]
    for idx1, idx2 in index_pairs(len(edge_sets)):
        print names[idx1], names[idx2]
        jaccard_index(edge_sets[idx1], edge_sets[idx2])


def all_pairs(path_sets, names, colors):
    edge_sets = [set(get_pairs_of_edges(paths)) for paths in path_sets]
    for idx1, idx2 in index_pairs(len(edge_sets)):
        print names[idx1], names[idx2]
        jaccard_index(edge_sets[idx1], edge_sets[idx2])


def all_identical(path_sets, names, colors):
    path_sets = [set(map(tuple, paths)) for paths in path_sets]
    for idx1, idx2 in index_pairs(len(path_sets)):
        ident = len(path_sets[idx1].intersection(path_sets[idx2]))
        print names[idx1], names[idx2]
        print "Identical:", ident, \
            "S1 percent:", ident * 1.0 / len(path_sets[idx1]), "S2 percent:", ident * 1.0 / len(path_sets[idx2])

    fig, ax = plt.subplots(1, 2, figsize=(20, 8))

    identical_venn(path_sets[:3], names[:3], ax[0], colors[:3])
    identical_venn(path_sets[-3:], names[-3:], ax[1], colors[-3:])


def all_foldchange(sg, path_sets, names):
    #lfcs = []
    fig, ax = plt.subplots(3, 2, figsize=(20, 12))
    ax = [ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1], ax[2, 0]]
    for i in xrange(len(path_sets)):
        values = map(lambda path: math.log(foldchange(*get_path_count(sg, path)), 2), path_sets[i])
        #lfcs += values
        print "%s: '-Inf' values: %f, +Inf: %f" % (names[i], 1.*len([x for x in values if x <= -15])/len(values),
                                            1. * len([x for x in values if x >= 15])/len(values))
        #sns.distplot([x for x in values if -15 < x < 15], rug=False, kde=False, bins=100, ax=ax[i])
        sns.distplot(values, rug=False, kde=False, bins=100, ax=ax[i])

        ax[i].set(xlabel="log2(fold change)")
        sns.distplot([x for x in values if -15 < x < 15], rug=False, hist=False, kde=True,
                     label=names[i], ax=ax[len(names)])

    ax[len(names)].set(xlabel="log2(fold change)")
    sns.plt.legend()
    #df = pandas.DataFrame(data={'heuristic': h, 'lfc': lfcs})
    #return df


def one_foldchange(sg, path_sets, names, colors):
    #lfcs = []
    fig, ax = plt.subplots(1, figsize=(15, 9))
    for i in xrange(len(path_sets)):
        values = map(lambda path: math.log(foldchange(*get_path_count(sg, path)), 2), path_sets[i])
        #lfcs += values
        print "%s: '-Inf' values: %f, +Inf: %f" % (names[i], 1.*len([x for x in values if x <= -15])/len(values),
                                            1. * len([x for x in values if x >= 15])/len(values))

        sns.distplot([x for x in values if -12 < x < 12], rug=False, hist=False, kde=True,
                     label=names[i], ax=ax, color=colors[i])

    ax.set(xlabel="log2(fold change)")
    ax.set_xlim(-12,12)
    plt.legend()
    #df = pandas.DataFrame(data={'heuristic': h, 'lfc': lfcs})
    #return df
    return ax


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


def prepare_coverage_data(sg, path_sets, names, min_length=1):

    heuristic = []
    variance0 = []
    dispersion0 = []
    variance1 = []
    dispersion1 = []

    for paths, name in zip(path_sets, names):
        print "counting ", name
        for path in paths:
            if len(path) >= min_length:
                heuristic += [name]
                cov0, cov1 = sg.get_path_coverage(path)
                variance0.append(variance(cov0))
                variance1.append(variance(cov1))
                dispersion0.append(dispersion(cov0))
                dispersion1.append(dispersion(cov1))

    df = pandas.DataFrame(data={'heuristic': heuristic, "variance0": variance0, "variance1": variance1,
                                "dispersion0": dispersion0, "dispersion1": dispersion1})
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


def all_coverage_variance(df, names, colors, remove=0.05):
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


def all_blast(names, prefix, colors):
    qcovs, nhits = count_qcov(['blast_' + prefix + '_' + name + '.tsv' for name in names], 1000)
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
                for read in to_old_names[node].split('|'):
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




