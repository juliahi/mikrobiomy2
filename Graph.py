import numpy as np
from collections import defaultdict, Counter
import Queue as Q
import bisect

def compl(s):
        if s == '' or s == []: return ''
        if s == None: return None
        def c(x): 
            if x == 'A': return 'T' 
            elif x == 'C': return 'G' 
            elif x == 'G': return 'C'
            elif x == 'T': return 'A'
            return 'N'
        return ''.join(map(c, s ))

    
    
class Node:
    epsilon=1.
    def __init__(self, id, values, seq, k):
        self.id = id
        #self.cov = int(values[1])
        #self.Ocov = int(values[2])
        self.next_nodes = [] 
        self.prev_nodes = [] 
        self.nreads = [Node.epsilon, Node.epsilon] #defaultdict(0.)
        #self.reads = []
        self.n = int(values[0]) ##len(seq)
        
        self.selected = 0  #assembly number in which it was selected lastly
        
        self.seq = seq

        self.length = self.n + k - 1
        
        self.twin = None
        
    def add_reads(self, cond, reads):
            if type(reads) == list:
                self.nreads[cond] += len(reads)
                #self.reads += reads
            else:
                self.nreads[cond] += reads

    def foldChange(self):
        if self.nreads[1] == 0: return None
        return 1.*self.nreads[0] / self.nreads[1]
    
    
    def get_score(self):    # should return float!
        fc = self.foldChange()
        if fc > 1: fc = 1./fc
        return fc
        
        
    def get_fasta(self, k, s=0, e=0): #s,e - offsets from start/end, 0-based
        seq = ''
        n = self.n
        if s >= k-1: return self.seq[s-k+1:n-e]
        if e < n:  
            seq = self.seq[:n-e]
            e=n
        if e >= k-1: return compl(self.twin.get_fasta(k, e, s)[::-1]) + seq
    
        ### not whole sequence available!
        
        if s < n: seq_pocz = compl(self.twin.get_fasta(k, k-1, s)[::-1])
        else: seq_pocz = ''
        
        s = max(s, n) 
        
        #try without recursion
        for node in self.prev_nodes: 
            if node.length >= 2*(k-1) or node.n >= k-1-s  :
                missing = node.get_fasta(k, node.length-(k-1-s), e-n)
                return seq_pocz + missing + seq
        
        for node in self.next_nodes: 
            if node.length >= 2*(k-1) or node.n >= k-1-e :
                missing = node.get_fasta(k, s-n, node.length-(k-1-e))
                return seq_pocz + missing + seq

        #with recursion
        missing = ''
        for node in self.prev_nodes: 
                newmissing = node.seq[-(k-1-s):node.n-(e-n)]
                if len(newmissing) > len(missing):
                    missing = newmissing
                    
        e += len(missing)
        seq = missing + seq
        missing = ''
        for node in self.next_nodes:
                newmissing = node.get_fasta(k, s-n, k-1)  #only from twin seq
                if newmissing != None and len(newmissing) > len(missing):
                    missing = newmissing
        s += len(missing)
        seq_pocz += missing
        
        for node in self.prev_nodes: 
                missing = node.get_fasta(k, node.length-(k-1-s), e-n)
                if missing != None: return seq_pocz + missing + seq
        for node in self.next_nodes: 
                missing = node.get_fasta(k, s-n, node.length-(k-1-e))
                if missing != None: return seq_pocz + missing + seq
        return None
        
        

class FullNode(Node):

    def add_reads(self, cond,reads):
        #self.reads += [r.nr for r in reads]
        self.reads += reads
        self.nreads[cond] += len(reads)
        

class Arc:
    def __init__(self, values):
        self.start = int(values[0])
        self.end = int(values[1])
        self.mult = int(values[2])
        
        self.nreads = 0
        
class Read:
    def __init__(self, values):
        self.nr = int(values[0])
        self.node_start = int(values[1])
        self.read_start = int(values[2])
        

class Path:
    def __init__(self, first_node, assembly_id):
        self.nodes = [first_node]
        self.nodes_prev = []  #REVERSED list of nodes before self.nodes - extensions added to list end
        self.counts = first_node.nreads[:]
        self.length = first_node.length
        self.assembly_id = assembly_id
    def first(self):
        if self.nodes_prev == []: return self.nodes[0]
        return self.nodes_prev[-1]
    def last(self): return self.nodes[-1]
    def foldChange(self): return 1.0 * self.counts[0] / self.counts[1]
    def node_ids(self): 
        return [n.id for n in self.nodes_prev[::-1]] + [n.id for n in self.nodes]
    
    def get_nodes(self): 
        return self.nodes_prev[::-1] + self.nodes
    def name(self):
        return "%s_%d_%f"%(self.assembly_id, self.length, self.foldChange())
    
    
    
    def expand(self, borderFC):
        best = None
        prev = True
        fC = self.foldChange()
        bestfC = 1
        for node in self.first().prev_nodes:
            if node.selected == self.assembly_id: continue # loop! change if not including previously selected nodes! 
            newfC = (1.* self.counts[0] + node.nreads[0])/(1.* self.counts[1] + node.nreads[1])
            #print fC, borderFC, newfC
            if (fC > 1 and newfC > bestfC and newfC > borderFC) or (fC < 1 and newfC < bestfC and newfC > 1/borderFC):
                bestfC = newfC
                best = node
            
        for node in self.last().next_nodes:
            if node.selected == self.assembly_id: continue # loop! change if not including previously selected nodes! 
            newfC = (1.* self.counts[0] + node.nreads[0])/(1.* self.counts[1] + node.nreads[1])
            #print fC, borderFC, newfC
            if (fC > 1 and newfC > bestfC and newfC > borderFC) or (fC < 1 and newfC < bestfC and newfC > 1/borderFC):
                bestfC = newfC
                best = node
                prev = False
        if best == None: return False
    
        best.selected = self.assembly_id
        if not prev:
            self.nodes.append(best)
        else:
            self.nodes_prev.append(best)
        
        self.length += best.n
        self.counts = [self.counts[0] + best.nreads[0], self.counts[1] + best.nreads[1]]
        
        return True
    
    def select(self):
        pass
        #TODO - modify node counts or sth
        


def arr_id(node_id): 
        if node_id > 0:
            return 2*node_id-2
        else:
            return -2*node_id-1

class VelvetGraph:
        
    def node(self, node_id):
        return self.nodes[arr_id(node_id)]
        
        
    def add_reads(self,node_id, reads, sample_ids):
        node = self.node(node_id)
        
        
        ##for 2 conditions!
        s1 = sum([self.conds[x] for x in sample_ids])
        node.add_reads(0, len(sample_ids) - s1)
        node.add_reads(1, s1)
        
        ###if adding reads
        ##node.add_reads(0, [reads[i] for i in xrange(len(reads)) if conds[sample_id[i]] == 0])
        ##node.add_reads(1, [reads[i] for i in xrange(len(reads)) if conds[sample_id[i]] == 1])
        
        ###faster version, if not tracking reads
        ##cond_counts = Counter([self.conds[s] for s in sample_ids])
        ##for c,v in cond_counts.items():
                        ##node.add_reads(c,v)
    
    def add_arc(self, values):
        n1 = int(values[0])
        n2 = int(values[1])
        #a = Arc(values)
        #self.arcs.append(a)
        
        self.node(n1).next_nodes.append(self.node(n2))
        self.node(n2).prev_nodes.append(self.node(n1))
        
        self.node(-n2).next_nodes.append(self.node(-n1))
        self.node(-n1).prev_nodes.append(self.node(-n2))

    
    def read_full(self, filename): #TODO:zaktualizowac
        with open(filename) as f:
            self.n_nodes, self.n_seqs, self.k, _ = map(int, f.readline().split())
            self.nodes = [None] * (self.n_nodes+1)
            self.nodes[0] = FullNode([0,0,0,0], '', '')
            self.nodes[0].selected = -1
            for i in range(1,self.n_nodes+1):
                #assert f, "File corrupted"
                values = f.readline().strip().split()
                #assert ltype == "NODE", "file corrupted"
                self.nodes[i] = FullNode(values[1:], f.readline(), f.readline())
                
            self.arcs = []
            self.seqs = []
            self.reads = []

            while f:
                values = f.readline().strip().split()
                if values == []:
                    break
                ltype = values[0]
                values = values[1:]
                if ltype == "ARC":
                    self.add_arc(values)
                
                #if ltype == "SEQ":    
                    
                if ltype == "NR":
                    node_id = int(values[0])
                    n_reads = int(values[1])
                    reads = [Read(f.readline().strip().split()) for j in xrange(n_reads)]
                    self.reads += reads
                    self.add_reads(node_id, cond, reads)
                    
                        
            ##for a in self.arcs:
                ##a.nreads = self.count_reads(a.start, a.end)
    
    
    def __init__(self, filename, read_counts, conds): #2 conditions
        #read_counts: list of n. of reads from each sample
        #conds: list with 0 or 1 indicating which condition each sample represents
        
        self.conds = conds
        for i in xrange(1, len(read_counts)): read_counts[i] += read_counts[i-1]
        ###def get_cond(nr):
            ###i = len(read_counts)-1
            ###nr = int(nr)
            ###assert nr <= read_counts[i], 'more reads given to build graph than now'
            ###while i >= 0 and nr <= read_counts[i]: 
                ###i -= 1
            ###return conds[i+1]
        
        with open(filename) as f:
            self.n_nodes, self.n_seqs, self.k, _ = map(int, f.readline().split())
            self.nodes = np.ndarray((2*self.n_nodes,),dtype=np.object)
            #self.nodes[0] = Node(0, [0,0,0, 0, 0 ], '', self.k)   #dummy node
            #self.conodes = np.ndarray((self.n_nodes+1,),dtype=np.object)
            #self.conodes[0] = Node(0, [0,0,0, 0, 0 ], '',  self.k)   #dummy node
            
            for i in range(1,self.n_nodes+1):
                #assert f, "File corrupted"
                values = f.readline().rstrip('\n').split()
                #assert ltype == "NODE", "file corrupted"
                self.nodes[arr_id(i)] = Node(int(values[1]), values[2:], f.readline().rstrip('\n'), self.k)
                self.nodes[arr_id(-i)] = Node(-int(values[1]), values[2:], f.readline().rstrip('\n'), self.k)
                self.nodes[arr_id(i)].twin = self.nodes[arr_id(-i)]
                self.nodes[arr_id(-i)].twin = self.nodes[arr_id(i)]
                
                
            #self.arcs = []
            #self.seqs = []
            #self.reads = []
            
            while f:
                values = f.readline().rstrip().split()
                if values == []:
                    break
                ltype = values[0]
                values = values[1:]
                if ltype == "ARC":
                   self.add_arc(values)
                
                #if ltype == "SEQ":    
                    
                if ltype == "NR":
                    read_count = [0,0]
                    node_id = int(values[0])
                    n_reads = int(values[1])
                    
                    rids = [int(f.readline().rstrip().split()[0]) for j in xrange(n_reads)]
                    sample_ids = [bisect.bisect_left(read_counts, r) for r in rids]
                    self.add_reads(node_id, rids, sample_ids)

        #self.normalize()
            
    def count_reads(self, node1, node2):   #count reads going common to two nodes
        i,j,r=0,0,0
        reads1 = self.node(node1).reads
        reads2 = self.node(node2).reads
        while i<len(reads1) and j<len(reads2):
            if reads1[i] < reads2[j]: i += 1
            elif reads1[i] > reads2[j]: j += 1
            else:
                i += 1
                j += 1
                r += 1
        return r
    
    def get_fasta_ids(self, node_ids):
        seq = ''
        start = 0
        for node_id in node_ids:
            seq += self.node(node_id).get_fasta(self.k, start, 0) ##### + '.'
            start = self.k-1
        return seq
    
    def get_fasta(self, nodes):
        seq = ''
        start = 0
        for node in nodes:
            seq += node.get_fasta(self.k, start, 0) ##### + '.'
            start = self.k-1
        return seq
    
    
    def normalize(self):
        #TODO
        pass
      

    
    def start_queue(self):
        self.queue = Q.PriorityQueue(2*len(self.nodes)-2) #priority = foldChange or 1/foldChange 
        for n in self.nodes:
            self.queue.put((n.get_score(), n))
        
    def max_fc_node(self, assembly_id):
        #l = [(n.foldChange(), n) for n in self.nodes if n.selected == 0] + [(n.foldChange(), n) for n in self.conodes if n.selected == 0]
        node = None
        while not self.queue.empty():
            score, node = self.queue.get()
            if node.selected == 0: break
        if node != None: node.selected = assembly_id
        return node
    
    



