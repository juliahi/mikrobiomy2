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
            # assert start2 == 0, "%s %s %d %d %d %d" % (n1, n2, start1, end1, start2, end2)
            return n1, n2, start1, end1, start2, end2

    return None


def cons_pairs(input_list):
    for i in xrange(len(input_list)-1):
        yield input_list[i:i+2]


class SgaGraph:
    def __init__(self, conds=None, max_edges=128):
        self.nodes = {}
        self.counts = {}
        self.duplicates_dict = {}
        self.conds = conds
        self.graph = None  # to be initialized in load_graph
        self.filename = None
        self.max_edges = max_edges
        self.old_n_nodes = None

    def init_graph(self, filename):  # 1 or 2 conditions

        self.filename = filename
        # Leave only nodes with neighbours
        with open(filename) as f:
            # self.header = f.readline()
            # print self.header
            f.readline()

            noedges0, noedges1 = 0, 0
            while True:  # Nodes
                line = f.readline()
                if line[0] != 'V':
                    # print line
                    break
                node_id = line.split()[1]
                self.nodes[node_id] = (0, 0)
            print 'no. nodes in graph:\t', len(self.nodes)
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

            print 'no. edges in graph:\t', (noedges0 + noedges1)

        prevnodes = len(self.nodes)
        # original number of nodes in graph (with "lonely" disconnected nodes)
        self.old_n_nodes = prevnodes


        # if removing lonely vertices
        # self.nodes = {k: v for k, v in self.nodes.iteritems() if (v[0] > 0) or (v[1] > 0)}
        # print "Found %d disconnected nodes:\t%f of nodes" % (prevnodes - currnodes,
        #                                                      (prevnodes - currnodes) * 1. / prevnodes)
        # currnodes = len(self.nodes)
        # if not removing lonely vertices yet
        lonely = len([1 for v in self.nodes.itervalues() if (v[0] == 0) and (v[1] == 0)])
        print "Found %d disconnected nodes:\t%f of nodes" % (lonely,
                                                             lonely * 1. / prevnodes)
        print "Found %d >-< edges, \t%f of all edges" % (noedges1, noedges1 * 1. / (noedges0 + noedges1))


        sr = [v[0]+v[1] for k, v in self.nodes.iteritems() if self.is_super_repetitive(k)]
        print "Found %d super-repetitive nodes with %d super-repetitive edges" % (len(sr), sum(sr))

        # # super-repetitive for only one end of node:
        # srin = [v[0] for v in self.nodes.values() if v[0] > self.max_edges]
        # print "Found %d super-repetitive nodes with %d super-repetitive in-edges" % (len(srin), sum(srin))
        # srout = [v[1] for v in self.nodes.values() if v[1] > self.max_edges]
        # print "Found %d super-repetitive nodes with %d super-repetitive out-edges" % (len(srout), sum(srout))
        # print "%d super-repetitive nodes for both out- and in-edges" % (len([1 for v in self.nodes.values()
        #                                                                      if v[0] > self.max_edges
        #                                                                      and v[1] > self.max_edges]))

    def is_super_repetitive(self, node):
        v = self.nodes[node]
        return v[0] + v[1] > self.max_edges


    def load_graph(self, remove_super_repeats=True):
        # loads graph without: island nodes, >-< edges, super-repetitive edges if srrm=True
        self.graph = nx.DiGraph()

        with open(self.filename) as f:
            f.readline()

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
            print 'no. nodes loaded:\t', self.graph.number_of_nodes()

            dup_edges = 0
            if not remove_super_repeats:    # don't remove super-repetitive edges
                while line:  # Edges
                    edge = read_edge(line)
                    if edge is not None:
                        dup_edges += self.add_edge(edge)
                    line = f.readline()

            else:    # remove = don't load super-repetitive edges
                repetitive = {}
                while line:
                    edge = read_edge(line)
                    # if edge is not None and not self.is_super_repetitive(edge[0]) \
                    #         and not self.is_super_repetitive(edge[1]):
                    #     dup_edges += self.add_edge(edge)

                    # mimic SGA
                    if edge is not None:
                        if not self.is_super_repetitive(edge[0]) \
                                and not self.is_super_repetitive(edge[1]):
                            dup_edges += self.add_edge(edge)
                        else:
                            repetitive[edge[0]] = True
                            repetitive[edge[1]] = True
                    line = f.readline()

                rmedges = list(self.graph.edges(repetitive.keys()))
                print "Removing %d edges from %d super-repetitive nodes" % (len(rmedges), len(repetitive))
                self.graph.remove_edges_from(rmedges)



            print 'no. edges loaded:\t', self.graph.number_of_edges()
            print 'Found %d duplicated edges' % dup_edges



    def subgraph(self, nodes):
        newsg = SgaGraph(self.conds)
        newsg.filename = self.filename
        for node in nodes:
            newsg.counts[node] = self.counts[node]
            if node in self.duplicates_dict:
                newsg.duplicates_dict[node] = self.duplicates_dict[node]
        newsg.graph = self.graph.subgraph(nodes).copy()
        return newsg


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
        # returns 1 if edge was a duplicated edge

        # check if we can remove start2 and end1
        assert values[4] == 0, str(values)
        assert values[3]+1 == self.graph.nodes[values[0]]["length"], str(values)

        if self.graph.has_edge(values[0], values[1]):
            d = self.graph.get_edge_data(values[0], values[1])
            if d["end1"] - d["start1"] < values[3] - values[2]:
                self.graph[values[0], values[1]]['start1'] = values[2]
                self.graph[values[0], values[1]]['end1'] = values[3]
                # self.graph[values[0], values[1]]['start2'] = values[4]
                self.graph[values[0], values[1]]['end2'] = values[5]
                return 1
        if self.graph.has_node(values[0]) and self.graph.has_node(values[1]):
            self.graph.add_edge(values[0], values[1],
                                start1=values[2], end1=values[3],
                                # start2=values[4],
                                end2=values[5])
        return 0

    def add_dupl(self, name, count, reverse):
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
                        #print line
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
                            self.add_dupl(name, c, reverse)
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

    def node_degree_pairs(self):
        in_degrees = self.graph.in_degree()
        out_degrees = self.graph.out_degree()
        return [(v, out_degrees[k]) for k, v in in_degrees]

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

    def node_foldchange(self, node):
        return SgaGraph.foldchange(*self.counts[node])

    def remove_nodes(self, nodes):
        for v in nodes:
            try:
                del self.counts[v]
                # del self.nodes[v]
            except KeyError:
                pass
        self.graph.remove_nodes_from(nodes)

    def remove_edges(self, edges):
        self.graph.remove_edges_from(edges)

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

        self.remove_nodes(to_remove)
        return to_remove


    def in_length(self, node, node_len):
        # maximal length according to incoming edges
        val = 0
        for end in self.graph.in_edges(node, data="end2"):
            val = max(val, node_len - end - 1)
        return max

    def remove_deadends_by_length(self, minlength=200):
        # removes nodes shorter than minlength, with zero out-degree or in-degree
        # always use on simplified graph (after compress_simple_paths)
        to_remove = {}

        for node, indeg in self.graph.in_degree():
            if indeg == 0:
                to_remove[node] = 0
        for node, outdeg in self.graph.out_degree():
            if outdeg == 0:
                if node in to_remove: del to_remove[node]
                else: to_remove[node] = 1

        print "Found %d dead-ends" % len(to_remove)

        # use node length as dead-end length
        # for node, l in self.graph.nodes.data("length"):
        #     if node in to_remove:
        #         if l >= minlength:
        #             del to_remove[node]

        # use non-overlap length as dead-end length
        for node, l in self.graph.nodes.data("length"):
            if node in to_remove:
                if to_remove[node] == 1:    # out_deg=0
                    for n1, n2, end in self.graph.in_edges(node, data="end2"):
                        if l-end-1 >= minlength:
                            del to_remove[node]
                            break
                else:    # to_remove[node] == 0: #in_deg=0
                    for n1, n2, start in self.graph.out_edges(node, data="start1"):
                        if start >= minlength:
                            del to_remove[node]
                            break


        print "Remove %d short dead-ends" % len(to_remove)
        self.remove_nodes(to_remove.keys())
        return to_remove.keys()


    def compress_simple_paths(self):
        # print "Edges before compression:", self.number_of_edges()
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

        print "Compressing paths: Removing %d nodes, adding %d nodes (=number of simplified paths)" \
              % (sum([len(x) for x in paths]), len(paths))
        self._compress_paths_finish(paths)
        return paths

    def _compress_paths_finish(self, paths):
        for path in paths:
            if self.graph.has_node(path[-1]):
                new_nodes, new_edges = self.compress_path(path)

                # print "Adding %d nodes" % len(new_nodes)
                self.graph.add_nodes_from([new_nodes])
                self.graph.add_edges_from(new_edges)
                self.remove_nodes(path)

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
            edges_to_add.append((new_name, v2, {'start1': length + edge_d["start1"], 'end1': length + edge_d["end2"],
                                                # 'start2': edge_d["start2"],
                                                'end2': edge_d["end2"]}))
        for v1 in self.graph.predecessors(path[0]):
            for v, v2, edge_d in self.graph.edges(v1, data=True):
                if v2 == path[0]:
                    edges_to_add.append((v1, new_name, edge_d))

        last_dict = self.graph.node[path[-1]]
        length += last_dict["length"]
        return (new_name, {'length': length}), edges_to_add


    def remove_short_islands(self, minlength):
        node_list = []
        for node, deg in self.graph.degree():
            if deg == 0 and self.graph.node[node]["length"] < minlength:
                node_list.append(node)
        print "Removing %d short-island-nodes" % len(node_list)
        self.remove_nodes(node_list)

    def fix_super_repetitive(self, threshold=128):
        edge_list = []
        no_nodes_in = 0
        for node, deg in self.graph.in_degree():
            if deg >= threshold:
                edge_list += list(self.graph.in_edges(node))
                no_nodes_in += 1
        no_edges_in = len(edge_list)
        no_nodes_out = 0
        for node, deg in self.graph.out_degree():
            if deg >= threshold:
                edge_list += list(self.graph.out_edges(node))
                no_nodes_out += 1

        no_edges_out = len(edge_list) - no_edges_in
        print "Removing %d in-edges from %d nodes and %d out-edges from %d nodes" % (no_edges_in, no_nodes_in,
                                                                                     no_edges_out, no_nodes_out)
        self.remove_edges(edge_list)


    def get_nodes_sequence(self, nodes, graph_file=None):
        seqs = {}
        edges = {}
        for node in nodes:
            for read in node.split('|'):
                seqs[read] = None
            for n1, n2 in cons_pairs(node.split('|')):
                edges[(n1, n2)] = None

        if graph_file is None:
            graph_file = self.filename

        with open(graph_file) as f:
            f.readline()
            while True:  # Nodes
                line = f.readline()
                if line[0] != 'V':
                    break
                split = line.split()
                node_id = split[1]
                seq = split[2]
                if node_id in seqs: seqs[node_id] = seq
            print 'no. sequences loaded:', len(seqs)

            while line:  # Edges
                edge = read_edge(line)
                if edge is not None and (edge[0], edge[1]) in edges:
                    edges[(edge[0], edge[1])] = edge[2:]
                line = f.readline()
            print 'no. edges loaded:', len(edges)

        sequences = []
        for node in nodes:
            reads = node.split('|')
            seq = seqs[reads[0]]
            for read1, read2 in cons_pairs(reads):
                edge = edges[(read1, read2)]
                seq += seqs[read2][edge[3]+1:]

                assert seqs[read2][:edge[3]+1] == seqs[read1][edge[0]:], \
                    seqs[read2][:edge[3]+1] + ' ' + seqs[read1][edge[0]:]

            assert len(seq) == self.graph.node[node]["length"], "%d %d" % (len(seq), self.graph.node[node]["length"])
            sequences.append(seq)
        return sequences


    def save_simple_nodes(self, minlength, out_file, graph_file=None, remove=False):
        node_list = []
        for node, deg in self.graph.degree():
            if deg == 0 and self.graph.node[node]["length"] >= minlength:
                node_list.append(node)
                if len(node_list) % 10000 == 0:
                    print len(node_list)
                    # break
        print "Saving %d sequences" % len(node_list)

        if node_list == []:
            return []

        seqs = self.get_nodes_sequence(node_list, graph_file)

        counts = 0
        length = 0
        # writing to FASTA
        with open(out_file, 'w+') as f:
            for node, seq in zip(node_list, seqs):
                counts += self.counts[node][0] + self.counts[node][1]
                length += len(seq)
                f.write('>%s %d %d %f\n%s\n' % (node, self.counts[node][0], self.counts[node][1],
                                                self.node_foldchange(node), seq))
        print "Saved %d sequences with %d counts of total length %d" % \
              (len(node_list), counts, length)

        # removing nodes
        if remove:
            self.remove_nodes(node_list)
        return node_list, seqs


def get_largest_components(sg, ncomps):
    components = sorted(list(nx.weakly_connected_components(sg.graph)), key=lambda x: len(x), reverse=True)[:ncomps]
    components = [item for sublist in components for item in sublist]
    return sg.subgraph(components)



