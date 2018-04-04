# import Queue as Q
# import bisect
import math
import networkx as nx


def compl(s):
    if s == '' or s == []: return ''
    if s is None: return None

    def c(x):
        if x == 'A':
            return 'T'
        elif x == 'C':
            return 'G'
        elif x == 'G':
            return 'C'
        elif x == 'T':
            return 'A'
        return 'N'

    return ''.join(map(c, s))


# class Edge:
# def __init__(self, values, node1, node2):
# #self.name1 = values[0]
# self.start1 = int(values[2])
# self.end1 = int(values[3])
# #self.len1 = int(values[4])
# #self.name2 = values[1]
# self.start2 = int(values[5])
# self.end2 = int(values[6])
# #self.len2 = int(values[7])

# self.orient = int(values[8]) #(0=same, 1=reverse)
# #self.diffs = int(values[9])

# def whattype(self):
# if self.node1.length == abs(self.end1 - self.start1)+1: return 'incl1'
# if self.node2.length == abs(self.end2 - self.start2)+1: return 'incl2'
# if self.orient == 0:
# if self.start2 < self.end2: return '>->'
# else: return '<-<'
# else:
# if self.start2 < self.end2: return '>-<'
# else: return '<->'

def read_edge(line):
    line = line.split()[1:]
    if int(line[8]) == 0:   # not changing direction
            n1 = line[0]
            n2 = line[1]
            start1 = int(line[2])
            end1 = int(line[3])
            # #self.len1 = int(values[4])
            start2 = int(line[5])
            end2 = int(line[6])
            # #self.len2 = int(values[7])
            if start1 == 0:  # flip nodes
                n1, n2 = n2, n1
                start1, start2 = start2, start1
                end1, end2 = end2, end1
            assert start2 == 0, "%s %s %d %d %d %d" % (n1, n2, start1, end1, start2, end2)
            return n1, n2, start1, end1, start2, end2

    return None



class SgaGraph:

    def __init__(self, filename, conds=None):  # 1 or 2 conditions
        self.conds = conds
        self.nodes = {}
        self.counts = None
        self.duplicates_dict = {}
        self.graph = None  # to be initialized in load_graph

        # Leave only nodes with neighbours
        with open(filename) as f:
            self.header = f.readline()
            print self.header

            noedges0, noedges1 = 0, 0
            while True:  # Nodes
                line = f.readline()
                if line[0] != 'V':
                    # print line
                    break
                node_id = line.split()[1]
                self.nodes[node_id] = (0, 0)
            print 'no. nodes loaded:', len(self.nodes)
            while line:  # Edges
                edge = read_edge(line)
                if edge is not None:
                    n1 = self.nodes[edge[0]]
                    self.nodes[edge[0]] = (n1[0], n1[1] + 1)
                    n2 = self.nodes[edge[1]]
                    self.nodes[edge[1]] = (n2[0] + 1, n2[1])
                    noedges0 += 1
                else:
                    noedges1 += 1
                line = f.readline()

            print 'no. edges loaded:', (noedges0 + noedges1)

        prevnodes = len(self.nodes)
        # original number of nodes in graph (with "lonely" disconnected nodes)
        self.old_n_nodes = prevnodes

        self.nodes = {k: v for k, v in self.nodes.iteritems() if (v[0] > 0) or (v[1] > 0)}
        currnodes = len(self.nodes)

        print "Removed %d not connected nodes: %f of nodes" % (prevnodes - currnodes,
                                                               (prevnodes - currnodes) * 1. / prevnodes)
        print "Found %d >-< edges, %f of all edges" % (noedges1, noedges1 * 1. / (noedges0 + noedges1))



    def load_graph(self, filename):
        self.graph = nx.DiGraph()
        with open(filename) as f:
            self.header = f.readline()
            while True:  # Nodes
                line = f.readline()
                if line[0] != 'V':
                    # print line
                    break
                line = line.split()
                node_id = line[1]
                # # tags = line[3:]
                if node_id in self.nodes:
                    self.add_node(node_id, line[2])
            print 'no. nodes loaded:', self.graph.number_of_nodes()
            while line:  # Edges
                edge = read_edge(line)
                if edge is not None:
                    self.add_edge(edge)
                line = f.readline()

            print 'no. edges loaded:', self.graph.number_of_edges()

    def whatcond(self, node_id):
        cond = node_id.split(':')[-1]
        if cond[-2] == '/': cond = cond[:-2]
        # TODO: kontrola bledu?
        return self.conds[cond]

    def get_node(self, node_id):
        try:
            return self.graph.node[node_id]
        except KeyError:
            try:
                return self.graph.node[self.duplicates_dict[node_id]]
            except KeyError:
                return None  # node removed because had no connected nodes

    def add_node(self, node_id, seq):
        cond = self.whatcond(node_id)
        nreads = [0, 0]
        nreads[cond] = 1
        self.graph.add_node(node_id, length=len(seq))  # , nreads=nreads)  # , seq=seq)
        # assert self.graph.node[node_id]["length"] is not None, \
        #    str(self.graph.node[node_id]) + 'node_id + ' ' + seq + len(seq)'

    def add_edge(self, values):
        if self.graph.has_node(values[0]) and self.graph.has_node(values[1]):
            self.graph.add_edge(values[0], values[1],
                                start1=int(values[2]), end1=int(values[3]),
                                start2=int(values[4]), end2=int(values[5]))

    def add_dupl2(self, name, count, reverse):
        try:
            if reverse:
                self.counts[name][1 - self.whatcond(name)] += count
            else:
                self.counts[name][self.whatcond(name)] += count
        except KeyError:
            try:
                if reverse:
                    self.counts[self.duplicates_dict[name]][1 - self.whatcond(name)] += count
                else:
                    self.counts[self.duplicates_dict[name]][self.whatcond(name)] += count
            except KeyError:
                pass  # node removed because had no connected nodes

    def add_duplicates_asqg(self, filename):
        # filename = asqg file of overlap btw duplicates 
        # create dictionary of removed reads: existing reads
        with open(filename) as f:
            f.readline()
            line = ''
            countdups = 0
            try:
                while True:  # Nodes
                    line = f.readline()
                    if line[0] != 'V':
                        print line
                        break
                while True:  # Edges
                    if line and line[0] == 'E':
                        line = line.split()[1:]
                        name1 = line[0].split(',')[0]  # remove seqrank=X from name of duplicated seq
                        # if name1 in self.nodes:
                        if (line[2] == '0' and int(line[3]) == int(line[4]) - 1) or (
                                line[5] == '0' and int(line[6]) == int(line[7]) - 1):
                            if self.whatcond(name1) != self.whatcond(line[1]):
                                countdups += 1
                                self.duplicates_dict[name1] = line[1]
                    else:
                        break
                    line = f.readline()
            except Exception as e:
                print e, 'line:', line, "."
                pass
            print "Found %d duplicates" % countdups

    def add_duplicates_fasta_todict(self, filename, reverse=False):
        # reverse - add to the second condition because duplicate is from another condition
        # filename = fasta file after rmdup (not duplicated) -- contains counts
        if not self.counts:
            self.counts = {}
            for k in self.nodes:
                self.counts[k] = [0, 0]
                self.counts[k][self.whatcond(k)] = 1

        with open(filename) as f:
            added = 0
            addedc = 0
            line = ''
            try:
                while True:
                    line = f.readline()
                    if line and line[0] == '>':
                        _, name, number = line.split()
                        c = int(number.split('=')[1]) - 1
                        if c > 0:
                            added += 1
                            addedc += c
                            self.add_dupl2(name, c, reverse)
                        f.readline()
                    else:
                        break
            except Exception as e:
                print e.message, 'line:', line
                pass
            print 'added %d reads to %d sequences' % (addedc, added)

    # def reset_counts(self):
    #     for node_id, node in self.graph.nodes.data():
    #         node["nreads"] = [0, 0]
    #         cond = self.whatcond(node_id)
    #         node["nreads"][cond] = 1

    def simple_cycles(self):
        return len(list(nx.simple_cycles(self.graph)))

    def find_cycle(self):
        return list(nx.find_cycle(self.graph))

    def node_degrees(self):
        return [v for k, v in self.graph.degree()]

    def number_of_reads(self):
        return sum([sum(v) for v in self.counts.itervalues()])

    def number_of_nodes(self):
        return self.graph.number_of_nodes()

    def number_of_edges(self):
        return self.graph.number_of_edges()

    def get_lengths(self):
        return [x[1] for x in self.graph.nodes.data("length")]

    @staticmethod
    def foldchange(n1, n2):
        if n2 == 0: return float("inf")
        return 1. * n1 / n2

    """ NREADS using nodes attributes """
    # def node_nreads_sums(self):
    #     return [sum(v[1]) for v in self.graph.nodes.data("nreads")]
    #
    # def node_nreads_cond(self, cond):
    #     return [v[1][cond] for v in self.graph.nodes.data("nreads")]
    #
    # def node_nreads(self):
    #     return [v[1] for v in self.graph.nodes.data("nreads")]

    # def foldchanges(self):
    #     return [SgaGraph.foldchange(*v[1]) for v in self.graph.nodes.data("nreads")]
    #
    # def log2foldchanges(self):
    #     return [math.log(SgaGraph.foldchange(*v[1]), 2) if v[1][0] != 0 else float("-inf")
    #                                       for v in self.graph.nodes.data("nreads")]

    """ NREADS using counts dictionary  """

    def node_nreads_sums(self):
        return [sum(v) for v in self.counts.itervalues()]

    def node_nreads_cond(self, cond):
        return [v[cond] for v in self.counts.itervalues()]

    def node_nreads(self):
        return self.counts.values()

    def foldchanges(self):
        return [SgaGraph.foldchange(*v) for v in self.counts.itervalues()]

    def log2foldchanges(self):
        return [math.log(SgaGraph.foldchange(*v), 2) if v[0] != 0 else float("-inf") for v in self.counts.itervalues()]


    # def traverse(self, nodename, visited):
    # q = Q.Queue() # queue with nodes
    # # if not visited[nodename]:
    # q.put(self.nodes[nodename])
    # visited[nodename] = True
    # size = 1
    # while not q.empty():
    # node = q.get()
    # for neighbour in node.neighbours():
    # if not visited[neighbour.id]:
    # visited[neighbour.id] = True
    # q.put(neighbour)
    # size += 1

    # return size

    def connected_components(self):
        cc = list(nx.weakly_connected_components(self.graph))
        return len(cc), [len(s) for s in cc]

    def remove_tips(self):
        # TODO remove tips that are first in pair?
        to_remove = []
        for node, neighbours in self.graph.adjacency():
            tmp_remove = []
            for v in neighbours:
                if self.graph.degree(v) == 1:
                    tmp_remove.append(v)
            if tmp_remove and len(tmp_remove) < len(neighbours):
                to_remove += tmp_remove
        return to_remove
        # self.graph.remove_nodes_from(to_remove)

    def remove_tips_finish(self, to_remove):
        self.graph.remove_nodes_from(to_remove)
        for v in to_remove:
            del self.counts[v]
            del self.nodes[v]

    def compress_simple_paths(self):
        print "Edges before compression:", self.number_of_edges()

        paths = []
        for node, neighbours in self.graph.adjacency():
            if len(neighbours) == 1 and self.graph.in_degree(node) != 1:
                # build path
                path = [node]
                tmp_node = neighbours.keys()[0]

                while True:
                    if self.graph.in_degree(tmp_node) != 1:
                        break
                    path.append(tmp_node)
                    if self.graph.out_degree(tmp_node) != 1:
                        break
                    tmp_node = list(self.graph.neighbors(tmp_node))[0]
                if len(path) > 1:
                    paths.append(path)

        print "Removing %d nodes, adding %d nodes (=number of simplified paths)" \
              % (sum([len(x) for x in paths]), len(paths))
        return paths

    def compress_paths_finish(self, paths):
        for i in xrange(0, len(paths), 2000):
            print i
            for path in paths[i:i+2000]:
                if self.graph.has_node(path[-1]):
                    new_nodes, new_edges = self.compress_path(path)
                    for v in path:
                        del self.counts[v]
                        #del self.nodes[v]

                    #print "Adding %d nodes" % len(new_nodes)
                    self.graph.add_nodes_from([new_nodes])
                    self.graph.add_edges_from(new_edges)
                    self.graph.remove_nodes_from(path)

    def compress_path(self, path):
        new_name = '|'.join(path)

        counts = [0, 0]
        for node in path:
            c = self.counts[node]
            counts[0] += c[0]
            counts[1] += c[1]
        self.counts[new_name] = counts

        # TODO retrieve sequence
        # node_dict = self.graph[path[0]]
        length = 0  # node_dict["length"]

        for n1 in path[:-1]:
            _, n2, edgedict = list(self.graph.edges(n1, data=True))[0]
            length += edgedict["start1"]

        edges_to_add = []
        for v1, v2, edge_d in self.graph.edges(path[-1], data=True):
            edges_to_add.append((new_name, v2, {'start1': edge_d["start1"], 'end1': edge_d["end2"],
                                                'start2': length + edge_d["start2"], 'end2': length + edge_d["end2"]}))
        for v1 in self.graph.predecessors(path[0]):
            for v, v2, edge_d in self.graph.edges(v1, data=True):
                if v2 == path[0]:
                    edges_to_add.append((v1, new_name, edge_d))

        last_dict = self.graph.node[path[-1]]
        length += last_dict["length"]
        return (new_name, {'length': length}), edges_to_add

