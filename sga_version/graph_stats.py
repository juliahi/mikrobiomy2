# graph analysis for SGAFilter graph

import numpy as np
import pandas
import networkx as nx
# import SGAFilter


def short_summary(sg):
    print "Number of nodes", sg.number_of_nodes()
    print "Number of edges", sg.number_of_edges()
    print "Number of reads", sg.number_of_reads()


def simple_cycles(sg):
    return len(list(nx.simple_cycles(sg.graph)))


def find_cycle(sg):
    return list(nx.find_cycle(sg.graph))


def analyze_counts(sg):
    counts = sg.node_nreads()
    n_nodes = len(counts)
    count_sums = [sum(v) for v in counts]
    n_reads = sum(count_sums)

    print "largest read counts:", sorted(count_sums)[-10:]
    print "Condition 0: %d, condition1: %d, all: %d" % (
        sum([x[0] for x in counts]), sum([x[1] for x in counts]), n_reads)
    # print 'duplicates', n_reads - n_nodes

    print "% read counts > 50:", len([1 for x in count_sums if x > 50]) * 1.0 / n_nodes

    countsnon0 = [n[0]+n[1] for n in counts if n[0] > 0 and n[1] > 0]
    print "nodes with coverage on both strands:%d, %f, having %f of reads" \
          % (len(countsnon0), 1.*len(countsnon0) / n_nodes, 1.*sum(countsnon0) / n_reads)

    return count_sums


def analyze_count_pairs(sg):
    pairs = sg.node_nreads()
    single_read_nodes = pairs.count([1, 0]) + pairs.count([0, 1])
    print "Number of single read nodes: %d, fraction of nodes: %f" % \
          (single_read_nodes, single_read_nodes*1.0/len(pairs))
    df = pandas.DataFrame({'c0': [min(10, v[0]) for v in pairs],
                           'c1': [min(10, v[1]) for v in pairs],
                           })
    return df


def analyze_degrees(sg):
    # degrees = sg.node_degrees()
    # degrees0 = [x for x in degrees if x > 0]
    n_nodes = sg.number_of_nodes()
    pairs = sg.node_degree_pairs()
    sumdegrees = [x[0] + x[1] for x in pairs]
    sumdegrees0 = [x for x in sumdegrees if x > 0]
    num_monobranch = len([1 for x in pairs if ((x[0] > 1) and (x[1] <= 1)) or ((x[1] > 1) and (x[0] <= 1))])
    num_dibranch = len([1 for x in pairs if x[0] > 1 and x[1] > 1])
    num_simple = len([1 for x in pairs if x[0] == 1 or x[1] == 1])

    print "\nVertices: %d\tEdges: %d \tIslands: %d \tTips: %d \tMonobranch: %d \tDibranch: %d \tSimple: %d" % \
          (n_nodes, sum(sumdegrees0)/2, n_nodes-len(sumdegrees0), sumdegrees0.count(1),
           num_monobranch, num_dibranch, num_simple)
    print "Node degrees:"
    print '% zeros:', 1. * (n_nodes - len(sumdegrees0)) / n_nodes
    print '% ones:', 1. * sumdegrees0.count(1) / n_nodes
    print '% twos:', 1. * sumdegrees0.count(2) / n_nodes
    print 'maximal values:', sorted(sumdegrees0)[-10:]
    print '% large values (>= 10):', 1. * (len([1 for x in sumdegrees0 if x >= 10])) / n_nodes
    print 'mean value:', np.mean(sumdegrees)
    print 'median value:', np.median(sumdegrees)
    print ''
    return sumdegrees


def connected_components(sg):
    cc = list(nx.weakly_connected_components(sg.graph))
    return len(cc), [len(s) for s in cc]


def analyze_connected_components(sg):
    no, sizes = connected_components(sg)
    n_nodes = sg.number_of_nodes()
    print "Number of components:", no, "for nodes:", n_nodes
    c1, c2 = sizes.count(1), sizes.count(2)
    print "Single-node components: %d, fraction of components: %f, fraction of nodes: %f" \
          % (c1, 1. * c1 / no, c1 * 1. / n_nodes)
    print "Two-node components: %d, fraction of components: %f, fraction of nodes: %f" \
          % (c2, 1. * c2 / no, c2 * 2. / n_nodes)
    larges = [x for x in sizes if x >= 100]
    print "Large components (>=100 nodes): %d, fraction of components: %f, fraction of nodes: %f" \
          % (len(larges), 1. * len(larges) / no, sum(larges) * 1. / n_nodes)
    print ''

    return sizes


def analyze_foldchanges(sg):
    lfcs = sg.log2foldchanges()
    lfcs2 = [min(max(x, -100), 100) for x in lfcs]

    return lfcs2


def analyze_lengths(sg):
    lengths = sg.get_lengths()
    print 'maximal lengths:', sorted(lengths)[-10:]
    print 'minimal lengths:', sorted(lengths)[:10]
    print '% large lengths (>= 200):', 1. * (len([1 for x in lengths if x >= 200])) / len(lengths)
    print 'mean value:', np.mean(lengths)
    print 'median value:', np.median(lengths)
    print ''

    return lengths


def get_path_count(graph, path, new_names=False):
        if new_names:
            path = [graph.new_names[x] for x in path]
        counts_tmp0 = 0
        counts_tmp1 = 0
        for node in path:
            counts_tmp0 += graph.counts[node][0]
            counts_tmp1 += graph.counts[node][1]
        return counts_tmp0, counts_tmp1


def get_path_coverage(graph, path, new_names=False):
    if new_names:
        path = [graph.new_names[x] for x in path]
    # coverage of path from graph = reads will be extended to whole node

    cov0 = [self.counts[path[0]][0]] * self.graph.node[path[0]]["length"]
    cov1 = [self.counts[path[0]][1]] * self.graph.node[path[0]]["length"]

    for node1, node2 in cons_pairs(path):
        if (node1, node2) in self.graph.edges:
                edge = self.graph.edges[node1, node2]
                c0, c1 = self.counts[node2]

                l2 = self.graph.node[node2]["length"]
                for i in xrange(edge['end2']):
                    cov0[-1-i] += c0
                    cov1[-1-i] += c1

                cov0 += [c0] * (l2 - edge['end2'])
                cov1 += [c1] * (l2 - edge['end2'])

        else:
                c0, c1 = self.counts[node2]
                l2 = self.graph.node[node2]["length"]
                cov0 += [c0] * l2
                cov1 += [c1] * l2
    return cov0, cov1


def get_largest_components(sg, ncomps):
    components = sorted(list(nx.weakly_connected_components(sg.graph)), key=lambda x: len(x), reverse=True)[:ncomps]
    components = [item for sublist in components for item in sublist]
    return sg.subgraph(components)


def get_random_components(sg, ncomps):
    components = list(nx.weakly_connected_components(sg.graph))[:ncomps]
    components = [item for sublist in components for item in sublist]
    return sg.subgraph(components)
