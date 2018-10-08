
from common import *
from SGAFilter import get_path_count, PATH_SEP
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import numpy


import common


def jaccard_index(set1, set2):
    wspolne = len(set1.intersection(set2))
    print "Set1\set2=", len(set1) - wspolne, '\t', "Set1 n Set2=", wspolne, '\t', "Set2\set1=", len(set2) - wspolne
    return 1. * wspolne / len(set1.union(set2))


def get_paths_from_file(filename, sg, min_len = 0):
    paths = []
    seqs = []
    for name, seq in iter_fasta(filename):
        #print name.split()[1][8:].split(';')[:-1]
        path = name.split()[1][8:].split(';')[:-1]
        path = [sg.new_names[x.strip('+-')] for x in path]
        if len(seq) > min_len:
            paths.append(path)
            seqs.append(seq)

    return paths, seqs


def get_edges_from_file(filename):
    all_edges = []
    for line in open(filename):
        if line.startswith(">"):
            name = line.split()[0][1:].replace('&', '|')
            edges = set(cons_pairs(name.split('|')))
            all_edges += edges
    return set(all_edges)


def get_edges_from_list(paths):
    all_edges = []
    for names in paths:
        name = '|'.join(names)
        edges = set(cons_pairs(name.split('|')))
        all_edges += edges
    return set(all_edges)


def compare_used_edges(file1, file2):
    if type(file1) is str:
        edges1 = get_edges_from_file(file1)
    else:
        edges1 = get_edges_from_list(file1)
    if type(file2) is str:
        edges2 = get_edges_from_file(file2)
    else:
        edges2 = get_edges_from_list(file2)

    index = jaccard_index(edges1, edges2)
    print "#edges1=%d\t#edges2=%d\tJaccard index=%f" % (len(edges1), len(edges2), index)


def get_edges_from_list_new(paths):
    all_edges = []
    for names in paths:
        edges = set(cons_pairs(names))
        all_edges += edges
    return set(all_edges)


def compare_used_edges_new(file1, file2):
    edges1 = get_edges_from_list_new(file1)
    edges2 = get_edges_from_list_new(file2)

    index = jaccard_index(edges1, edges2)
    print "#edges1=%d\t#edges2=%d\tJaccard index=%f" % (len(edges1), len(edges2), index)


def compare_used_vertices(paths1, paths2):
    nodes1 = set([node for path in paths1 for node in path ])
    nodes2 = set([node for path in paths2 for node in path ])

    index = jaccard_index(nodes1, nodes2)
    print "#vertices1=%d\t#vertices2=%d\tJaccard index=%f" % (len(nodes1), len(nodes2), index)


def vertices_venn(paths, names):
    v1 = [node for path in paths[0] for node in path]
    v2 = [node for path in paths[1] for node in path]
    v3 = [node for path in paths[2] for node in path]

    venn3([set(v1), set(v2), set(v3)],  set_labels=names)
    plt.title('Comparison of used vertices\n')


def get_pairs_of_edges(paths):
    all_pairs = []
    for path in paths:
        if len(path) >= 3:
            all_pairs += [(path[i], path[i+1], path[i+2]) for i in xrange(len(path)-2)]
    return set(all_pairs)


def compare_pairs_of_edges(paths1, paths2):
    edges1 = get_pairs_of_edges(paths1)
    edges2 = get_pairs_of_edges(paths2)

    index = jaccard_index(edges1, edges2)
    print "#pairs1=%d\t#pairs2=%d\tJaccard index=%f" % (len(edges1), len(edges2), index)


def compare_identical(paths1, paths2):
    s1 = set(map(tuple, paths1))
    s2 = set(map(tuple, paths2))
    ident = len(s1.intersection(s2))

    print "Identical:", ident, "S1 percent:", ident*1.0/len(s1), "S2 percent:", ident*1.0/len(s2)


def identical_venn(paths, names):
    venn3([set(paths[0]), set(paths[1]), set(paths[2])],
          set_labels=names)
    plt.title('Comparison of returned paths\n')


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


def stats(graph, paths, seqs):
    print "No. paths:", len(paths), "No. vertices:", sum([len(path) for path in paths]), \
        "Length (bp): ", sum(map(len, seqs)), \
        "Counts:", sum(map(lambda x:sum(get_path_count(graph, x)), paths))


def get_paths_fc(sg, paths, single=True):
    #single - include single-node paths
    if single:
        fcs = [map(sg.node_foldchange, nodes) for nodes in paths]
        lens = [map(lambda x:sg.graph.nodes[x]["length"], nodes) for nodes in paths]
    else:

        fcs = [map(sg.node_foldchange, nodes) for nodes in paths if len(nodes) > 1]
        lens = [map(lambda x:sg.graph.nodes[x]["length"], nodes) for nodes in paths if len(nodes) > 1]
    return fcs, lens


def fc_variance(fc, total_lfc):
    return sum([(math.log(x, 2)-total_lfc)**2 for x in fc])/len(fc)
    #return numpy.var([math.log(x, 2) for x in fc])
    #return numpy.var(fc)


def fc_weighted_variance(fc, weights):
    s = 1.*sum(weights)
    m = numpy.mean([math.log(x, 2)*y/s for x, y in zip(fc, weights)])
    return sum([y/s*((math.log(x, 2)*y/s - m)**2) for x, y in zip(fc, weights)])
    # return numpy.var([1.*x*y/s for x, y in zip(fc, weights)])


# def compare_lengths(file1, file2):
#     info1 = get_info_from_file(file1)
#     info2 = get_info_from_file(file2)
#
#     # TODO
#     # analyse lengths?

