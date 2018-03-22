import numpy as np
from collections import defaultdict, Counter
import Queue as Q
import bisect
import math
import networkx as nx


def compl(s):
    if s == '' or s == []: return ''
    if s == None: return None

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


# class Node:
##epsilon=.001
# def __init__(self, id, seq):
# self.id = id
# self.next_edges = []
# self.prev_edges = []
# self.nreads = [0, 0] #defaultdict(0.)
##self.n = int(values[0]) ##=len(seq)

##self.selected = 0  #assembly number in which it was selected lastly
#####self.seq = seq
# self.length = len(seq)


# def foldChange(self):
# if self.nreads[1] == 0: return float("inf")
# return 1.*self.nreads[0] / self.nreads[1]


# def get_score(self):    # should return float!
# fc = self.foldChange()
# if fc > 1: fc = 1./fc
# return fc

# def add_edge(self, edge, first):
# if first == 0:
# self.next_edges.append(edge)
# else:
# self.prev_edges.append(edge)

# def neighbours(self): #TODO: co jesli zmienia sie kierunek?
# return [e.node2 for e in self.next_edges] + [e.node1 for e in self.prev_edges]


# class Edge:
# def __init__(self, values, node1, node2):
##self.name1 = values[0]
# self.start1 = int(values[2])
# self.end1 = int(values[3])
##self.len1 = int(values[4])
##self.name2 = values[1]
# self.start2 = int(values[5])
# self.end2 = int(values[6])
##self.len2 = int(values[7])

# node1.add_edge(self, 0)
# node2.add_edge(self, 1)
# self.node1 = node1
# self.node2 = node2

# self.orient = int(values[8]) #(0=same, 1=reverse)
##self.diffs = int(values[9])

# def whattype(self):
# if self.node1.length == abs(self.end1 - self.start1)+1: return 'incl1'
# if self.node2.length == abs(self.end2 - self.start2)+1: return 'incl2'
# if self.orient == 0:
# if self.start2 < self.end2: return '>->'
# else: return '<-<'
# else:
# if self.start2 < self.end2: return '>-<'
# else: return '<->'

# def join_nodes(self):
# if len(self.node1.next_edges) == 1 and len(self.node2.prev_edges) == 1: #join
# newid = self.node1.id + '|' + self.node2.id
# if self.orient == 0:
##newseq = self.node1.seq += self.node2.seq[self.end2+1:]
# self.node1.next_edges = self.node2.next_edges
# x = self.start1
# for e in self.node1.next_edges:
# e.start1 += x
# e.end1 += x
# e.node1 = self.node1
# self.node1.id = newid
# self.node1.nreads = [c[0] + c[1] for c in zip(self.node1.nreads, self.node2.nreads)]
# self.node1.length += self.end2 - self.start2 + 1
# return True
# return False


class SgaGraph:

    def __init__(self, filename, conds=None):  # 1 or 2 conditions
        self.conds = conds
        self.nodes = {}
        self.duplicates_dict = {}
        self.graph = None   # to be initialized in load_graph

        ## Leave only nodes with neighbours
        with open(filename) as f:
            self.header = f.readline()
            print self.header

            noedges0, noedges1 = 0, 0
            while True:  # Nodes
                    line = f.readline()
                    if line[0] != 'V':
                        print line
                        break
                    node_id = line.split()[1]
                    self.nodes[node_id] = (0, 0)
            print 'no. nodes loaded:', len(self.nodes)
            while line:  # Edges
                    line = line.split()
                    n1 = self.nodes[line[1]]
                    self.nodes[line[1]] = (n1[0], n1[1] + 1)
                    n2 = self.nodes[line[2]]
                    self.nodes[line[2]] = (n2[0] + 1, n2[1])
                    if int(line[9]) == 0:
                        noedges0 += 1
                    else:
                        noedges1 += 1
                    line = f.readline()

            print 'no. edges loaded:', (noedges0 + noedges1)

        prevnodes = len(self.nodes)

        counter = 0
        for k, v in self.nodes.iteritems():
            if v[0] == 0 and v[1] == 0:
                print k
                counter += 1
            if counter > 100:
                break

        self.nodes = {k: v for k, v in self.nodes.iteritems() if (v[0] > 0)}
        currnodes = len(self.nodes)

        print "Removed %d not connected nodes: %f of nodes" % (
        prevnodes - currnodes, (prevnodes - currnodes) * 1. / prevnodes)
        print "Found %d >-< edges, %f of all edges" % (noedges1, noedges1 * 1. / (noedges0 + noedges1))

    def load_graph(self, filename):
        self.graph = nx.DiGraph()
        with open(filename) as f:
            self.header = f.readline()
            while True:  # Nodes
                line = f.readline()
                if line[0] != 'V':
                    print line
                    break
                line = line.split()
                node_id = line[1]
                ## tags = line[3:]
                if node_id in self.nodes:
                    self.add_node(node_id, line[2])
            print 'no. nodes loaded:', self.graph.number_of_nodes()
            while line:  # Edges
                line = line.split()[1:]
                if int(line[8]) == 0:
                    self.add_edge(line)
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
                return None

    def add_node(self, node_id, seq):
        cond = self.whatcond(node_id)
        nreads = [0, 0]
        nreads[cond] = 1
        self.graph.add_node(node_id, length=len(seq), nreads=nreads)  # , seq=seq)

    def add_edge(self, values):
        self.graph.add_edge(values[0], values[1],
                            start1=int(values[2]), end1=int(values[3]),
                            start2=int(values[5]), end2=int(values[6]))

    def add_dupl(self, name, count, reverse):
        node = self.get_node(name)
        if node:
            if reverse:
                node["nreads"][1 - self.whatcond(name)] += count
            else:
                node["nreads"][self.whatcond(name)] += count

    def add_duplicates_fasta(self, filename, reverse=False):
        # reverse - add to the second condition because duplicate is from another condition
        # filename = fasta file after rmdup (not duplicated) -- contains counts
        with open(filename) as f:
            try:
                added = 0
                addedc = 0
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
                print e, 'line:', line
                pass
            print 'added %d reads to %d sequences' % (addedc, added)

    def add_duplicates_asqg(self, filename):
        # filename = asqg file of overlap btw duplicates 
        # create dictionary of removed reads: existing reads
        with open(filename) as f:
            f.readline()
            try:
                while True:  # Nodes
                    line = f.readline()
                    if line[0] != 'V':
                        print line
                        break
                # lines = 0
                countdups = 0
                while True:  # Edges
                    if line and line[0] == 'E':
                        line = line.split()[1:]
                        # lines += 1
                        # if lines % 100000 == 0: print lines
                        name1 = line[0].split(',')[0]  # remove seqrank=X from name of duplicated seq
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

    def reset_counts(self):
        for node_id, node in self.graph.nodes.data():
            node["nreads"] = [0, 0]
            cond = self.whatcond(node_id)
            node["nreads"][cond] = 1

    def simple_cycles(self):
        return len(list(nx.simple_cycles(self.graph)))

    def find_cycle(self):
        return list(nx.find_cycle(self.graph))

    def node_degrees(self):
        return [v for k, v in self.graph.degree()]

    def node_nreads_sums(self):
        return [sum(v[1]) for v in self.graph.nodes.data("nreads")]

    def node_nreads_cond(self, cond):
        return [v[1][cond] for v in self.graph.nodes.data("nreads")]

    def node_nreads(self):
        return [v[1] for v in self.graph.nodes.data("nreads")]


    @staticmethod
    def foldchange(n1, n2):
        if n2 == 0: return float("inf")
        return 1.*n1/n2

    def foldchanges(self):
        return [SgaGraph.foldchange(*v[1]) for v in self.graph.nodes.data("nreads")]

    def log2foldchanges(self):
        return [math.log(SgaGraph.foldchange(*v[1]), 2) if v[1][0] != 0 else float("-inf") for v in self.graph.nodes.data("nreads")]

    # def edge_types(self):
    #     return [e.whattype() for e in self.edges]

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
        to_remove = []
        for node, neighbours in self.graph.adjacency():
            tmp_remove = []
            for v, _ in neighbours:
                if self.graph.degree(v)[1] == 1:
                    tmp_remove.append(v)
            if tmp_remove and len(tmp_remove) < len(neighbours):
                to_remove += tmp_remove
        return to_remove
        #self.graph.remove_nodes_from(to_remove)

    def compress_simple_paths(self):
        print "Edges before compression:", len(self.graph.number_of_edges())

        to_remove = []
        new_nodes = []
        new_edges = []
        for node, neighbours in self.graph.adjacency():
            if len(neighbours) == 1 and self.graph.in_degree(node)[1] != 1:
                # build path

                path = [node]
                tmp_node = neighbours.keys()[0]

                while True:
                    if self.graph.in_degree(tmp_node)[1] != 1:
                        break
                    path += tmp_node
                    tmp_node = self.graph.neighbours(tmp_node)[0]
                    if self.graph.out_degree(tmp_node)[1] != 1:
                        break

                if len(path) > 1:
                    tmp = self.compress_path(path)
                    new_nodes.append(tmp[0])
                    new_edges += tmp[1]
        to_remove += path

        print "Removing %d nodes" % len(to_remove)
        #return to_remove, new_nodes, new_edges
        self.graph.remove_nodes_from(to_remove)

        print "Adding %d nodes" % len(new_nodes)
        self.graph.add_nodes_from(new_nodes)
        self.graph.add_edges_from(new_edges)



    def compress_path(self, path):
        new_name = '|'.join(path)

        counts = [0, 0]
        for node, c in self.graph.nodes(path, data="counts"):
            counts[0] += c[0]
            counts[1] += c[1]

        # TODO retrieve sequence
        # node_dict = self.graph.nodes(path[0], data=True)[1]
        length = 0  # node_dict["length"]

        for n1 in path[:-1]:
            _, n2, edgedict = self.graph.edges(n1, data=True)[0]
            length += edgedict["start1"]

        edges_to_add = []
        for v1, v2, edge_d in self.graph.edges(path[-1], data=True):
            edges_to_add.append((new_name, v2, {'start1': edge_d["start1"], 'end1':edge_d["end2"],
                                'start2':length + edge_d["start2"], 'end2':length+edge_d["end2"]}))

        last_dict = self.graph.nodes(path[-1], data=True)[1]
        length += last_dict["length"]
        return (new_name, {'length': length+last_dict["length"]}), edges_to_add



    # def remove_edge(self, edge):
    #     try:
    #         edge.node1.next_edges.remove(edge)
    #     except:
    #         edge.node1.prev_edges.remove(edge)
    #
    #     try:
    #         edge.node2.prev_edges.remove(edge)
    #     except:
    #         edge.node2.next_edges.remove(edge)

