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
    epsilon=.1
    def __init__(self, id, seq):
        self.id = id
        self.next_nodes = [] 
        self.prev_nodes = [] 
        self.nreads = [Node.epsilon, Node.epsilon] #defaultdict(0.)
        #self.n = int(values[0]) ##=len(seq)
        
        self.selected = 0  #assembly number in which it was selected lastly
        ####self.seq = seq
        

    def foldChange(self):
        if self.nreads[1] == 0: return None
        return 1.*self.nreads[0] / self.nreads[1]
    
    
    def get_score(self):    # should return float!
        fc = self.foldChange()
        if fc > 1: fc = 1./fc
        return fc
        
    def add_edge(self, edge, first):
        if first == 0:
            self.next_nodes.append(edge)
        else:
            self.prev_nodes.append(edge)
        
    #def get_fasta(self, k, s=0, e=0): #s,e - offsets from start/end, 0-based
        
        #seq = ''
        #n = self.n
        #if n+k-1-e <= s: return ''
        #if s >= k-1: return self.seq[s-k+1:n-e]
        #if e < n:  
            #seq = self.seq[:n-e]
            #e=n
        #if e >= k-1: return compl(self.twin.get_fasta(k, e, s)[::-1]) + seq
    
        #### not whole sequence available!
        
        #if s < n: seq_pocz = compl(self.twin.get_fasta(k, k-1, s)[::-1])
        #else: seq_pocz = ''
        
        #s = max(s, n) 
        
        ##try without recursion
        #for node in self.prev_nodes: 
            #if node.length >= 2*(k-1) or node.n >= k-1-s  :
                #missing = node.get_fasta(k, node.length-(k-1-s), e-n)
                #return seq_pocz + missing + seq
        
        #for node in self.next_nodes: 
            #if node.length >= 2*(k-1) or node.n >= k-1-e :
                #missing = node.get_fasta(k, s-n, node.length-(k-1-e))
                #return seq_pocz + missing + seq

        ##with recursion
        #missing = ''
        #for node in self.prev_nodes: 
                #newmissing = node.seq[-(k-1-s):node.n-(e-n)]
                #if len(newmissing) > len(missing):
                    #missing = newmissing
                    
        #e += len(missing)
        #seq = missing + seq
        #missing = ''
        #for node in self.next_nodes:
                #newmissing = node.get_fasta(k, s-n, k-1)  #only from twin seq
                #if newmissing != None and len(newmissing) > len(missing):
                    #missing = newmissing
        #s += len(missing)
        #seq_pocz += missing
        
        #for node in self.prev_nodes: 
            #if node != self:
                #missing = node.get_fasta(k, node.length-(k-1-s), e-n)
                #if missing != None: return seq_pocz + missing + seq
        #for node in self.next_nodes: 
            #if node != self:
                #missing = node.get_fasta(k, s-n, node.length-(k-1-e))
                #if missing != None: return seq_pocz + missing + seq
        #return None



class Edge:
    def __init__(self, values):
        self.name1 = values[0]
        self.start1 = int(values[2])
        self.end1 = int(values[3])
        self.len1 = int(values[4])
        self.name2 = values[1]
        self.start2 = int(values[5])
        self.end2 = int(values[6])
        self.len2 = int(values[7])
        
        
        
        self.orient = int(values[8]) #(0=same, 1=reverse) 
        #self.diffs = int(values[9]) 
        
    def whattype(self):
        if self.len1 == abs(self.end1 - self.start1)+1: return 'incl1'
        if self.len2 == abs(self.end2 - self.start2)+1: return 'incl2'
        if self.orient == 0: 
            if self.start1 < self.end1: return '>->'
            else: return '<-<'
        else:
            if self.start1 < self.end1: return '>-<'
            else: return '<->'


#class Path:
    #def __init__(self, first_node, assembly_id):
        #self.nodes = [first_node]
        #self.nodes_prev = []  #REVERSED list of nodes before self.nodes - extensions added to list end
        #self.counts = first_node.nreads[:]
        #self.length = first_node.length
        #self.assembly_id = assembly_id
    #def first(self):
        #if self.nodes_prev == []: return self.nodes[0]
        #return self.nodes_prev[-1]
    #def last(self): return self.nodes[-1]
    #def foldChange(self): return 1.0 * self.counts[0] / self.counts[1]
    #def node_ids(self): 
        #return [n.id for n in self.nodes_prev[::-1]] + [n.id for n in self.nodes]
    
    #def get_nodes(self): 
        #return self.nodes_prev[::-1] + self.nodes
    #def name(self):
        #return "%s_%d_%f"%(self.assembly_id, self.length, self.foldChange())
    
    
    
    #def expand(self, borderFC):
        #best = None
        #prev = True
        #fC = self.foldChange()
        #bestfC = 1
        #for node in self.first().prev_nodes:
            #if node.selected == self.assembly_id: continue # loop! change if not including previously selected nodes! 
            #newfC = (1.* self.counts[0] + node.nreads[0])/(1.* self.counts[1] + node.nreads[1])
            ##print fC, borderFC, newfC
            #if (fC > 1 and newfC > bestfC and newfC > borderFC) or (fC < 1 and newfC < bestfC and newfC > 1/borderFC):
                #bestfC = newfC
                #best = node
            
        #for node in self.last().next_nodes:
            #if node.selected == self.assembly_id: continue # loop! change if not including previously selected nodes! 
            #newfC = (1.* self.counts[0] + node.nreads[0])/(1.* self.counts[1] + node.nreads[1])
            ##print fC, borderFC, newfC
            #if (fC > 1 and newfC > bestfC and newfC > borderFC) or (fC < 1 and newfC < bestfC and newfC > 1/borderFC):
                #bestfC = newfC
                #best = node
                #prev = False
        #if best == None: return False
    
        #best.selected = self.assembly_id
        #if not prev:
            #self.nodes.append(best)
        #else:
            #self.nodes_prev.append(best)
        
        #self.length += best.n
        #self.counts = [self.counts[0] + best.nreads[0], self.counts[1] + best.nreads[1]]
        
        #return True
    
    #def select(self):
        #pass
        ##TODO - modify node counts or sth
        


#def arr_id(node_id): 
        #if node_id > 0:
            #return 2*node_id-2
        #else:
            #return -2*node_id-1

class SgaGraph:
    def whatcond(self, node_id):
        cond = node_id.split(':')[-1]
        #TODO: kontrola bledu?
        return self.conds[cond]
    
    def node(self, node_id):
        return self.nodes[node_id]
    
    def add_edge(self, values):
        a = Edge(values)
        self.edges.append(a)  #?
        self.node(values[0]).add_edge(a, 0)
        self.node(values[1]).add_edge(a, 1)
        
    def add_dupl(self, values):
        name2 = values[1]
        self.node(values[0]).nreads[self.whatcond(name2)] += 1
        
    def add_node(self, node_id, seq):
        node = Node(node_id, seq)
        if self.conds:
            node.nreads[self.whatcond(node_id)] += 1
        self.nodes[node_id] = node
        
    def __init__(self, filename, conds=None): #1 or 2 conditions
        self.conds = conds
        self.nodes = {}
        self.edges = []
        with open(filename) as f:
            self.header = f.readline()
            print self.header
            
            try:
                while True: #Nodes
                    line = f.readline()
                    if line[0] != 'V':
                        print line
                        break
                    line = line.split()
                    node_id = line[1]
                    seq = line[2]
                    ##tags = line[3:]
                    self.add_node(node_id, seq)
                print 'no. nodes loaded:', len(self.nodes)
                while True: #Edges
                    line = line.split()[1:]
                    self.add_edge(line)
                    line = f.readline()
            except: pass
                    
    #def add_duplicates(self, filename): 
        #with open(filename) as f:
            #self.header = f.readline()
            #print self.header
            #cond = None
            #try:
                #while True: #Nodes
                    #line = f.readline()
                    #if line[0] != 'V':
                        #print line
                        #break
                #while True: #Edges
                    #line = line.split()[1:]
                    #self.add_dupl(line)
                    #line = f.readline()
            #except: pass
    
    
    def node_outdegrees(self):
        return [len(x.next_nodes)  for x in self.nodes.values()]
    def node_indegrees(self):
        return [len(x.prev_nodes)  for x in self.nodes.values()]
    def edge_types(self):
        return [e.whattype()  for e in self.edges]
    
    #def get_fasta_ids(self, node_ids):
        #seq = ''
        #start = 0
        #for node_id in node_ids:
            #seq += self.node(node_id).get_fasta(self.k, start, 0) ##### + '.'
            #start = self.k-1
        #return seq
    
    #def get_fasta(self, nodes):
        #seq = ''
        #start = 0
        #for node in nodes:
            #seq += node.get_fasta(self.k, start, 0) ##### + '.'
            #start = self.k-1
        #return seq
    
    






###########do testowania
def get_fasta(selfnode, k, s=0, e=0): #s,e - offsets from start/end, 0-based
        seq = ''
        self = selfnode
        n = self.n
        print self.id,s,e, n+k-1-s-e
        if n+k-1-e <= s: return ''
        if s >= k-1: return self.seq[s-k+1:n-e]
        if e < n:  
            seq = self.seq[:n-e]
            e=n
        if e >= k-1: return compl(get_fasta(self.twin, k, e, s)[::-1]) + seq
    
        ### not whole sequence available!
        
        if s < n: seq_pocz = compl(get_fasta(self.twin, k, k-1, s)[::-1])
        else: seq_pocz = ''
        
        s = max(s, n) 
        
        #try without recursion
        for node in self.prev_nodes: 
            if node.length >= 2*(k-1) or node.n >= k-1-s  :
                missing = get_fasta(node, k, node.length-(k-1-s), e-n)
                return seq_pocz + missing + seq
        
        for node in self.next_nodes: 
            if node.length >= 2*(k-1) or node.n >= k-1-e :
                missing = get_fasta(node, k, s-n, node.length-(k-1-e))
                return seq_pocz + missing + seq

        #with recursion
        
        print self.id, seq_pocz, seq, s, e
        
        missing = ''
        for node in self.prev_nodes: 
                newmissing = node.seq[-(k-1-s):node.n-(e-n)]
                if len(newmissing) > len(missing):
                    missing = newmissing
                    
        e += len(missing)
        seq = missing + seq
        missing = ''
        for node in self.next_nodes:
                newmissing = get_fasta(node, k, s-n, k-1)  #only from twin seq
                if newmissing != None and len(newmissing) > len(missing):
                    missing = newmissing
        s += len(missing)
        seq_pocz += missing
        
        print self.id, seq_pocz, seq
    
        for node in self.prev_nodes: 
            if node != self:
                print 'rec', node.id, node.length, node.length-(k-1-s), e-n
                #missing =  'N'*(node.length-e+n-node.length+(k-1-s)) #get_fasta(node, (k, node.length-(k-1-s), e-n)
                missing =  get_fasta_back(node, k, node.length-(k-1-s), e-n)
                if missing != None: return seq_pocz + missing + seq
        for node in self.next_nodes: 
            if node != self:
                print 'rec2', node.id, node.length, s-n, node.length-(k-1-e)
                #missing =  'N'*(node.length-s+n-node.length+(k-1-e))  ###get_fasta(node, (k, s-n, node.length-(k-1-e))
                missing =  get_fasta_forward(node, k, s-n, node.length-(k-1-e))
                if missing != None: return seq_pocz + missing + seq
        return None
        




def get_fasta_forward(selfnode, k, s=0, e=0): #s,e - offsets from start/end, 0-based
        seq = ''
        self=selfnode
        n = self.n
        if n+k-1-e <= s: return ''
        if s >= k-1: return self.seq[s-k+1:n-e]
        if e < n:  
            seq = self.seq[:n-e]
            e=n
        if e >= k-1: return compl(get_fasta(self.twin, k, e, s)[::-1]) + seq
    
        ### not whole sequence available!
        
        if s < n: seq_pocz = compl(get_fasta(self.twin, k, k-1, s)[::-1])
        else: seq_pocz = ''
        
        s = max(s, n) 
        for node in self.prev_nodes: 
            if node.length >= 2*(k-1) or node.n >= k-1-s  :
                missing = get_fasta_back(node, k, node.length-(k-1-s), e-n)
                return seq_pocz + missing + seq
        for node in self.next_nodes: 
            if node.length >= 2*(k-1) or node.n >= k-1-e :
                missing = get_fasta_forward(node, k, s-n, node.length-(k-1-e))
                return seq_pocz + missing + seq

        #with recursion
        missing = ''
        for node in self.next_nodes:
                newmissing = get_fasta_forward(node, k, s-n, k-1)  #only from twin seq
                if newmissing != None and len(newmissing) > len(missing):
                    missing = newmissing
        s += len(missing)
        seq_pocz += missing
        
        for node in self.next_nodes: 
            if node != self:
                missing = get_fasta_forward(node, k, s-n, node.length-(k-1-e))
                if missing != None: return seq_pocz + missing + seq
        return None
        
        

def get_fasta_back(selfnode, k, s=0, e=0): #s,e - offsets from start/end, 0-based
        seq = ''
        self = selfnode
        n = self.n
        print self.id,s,e, n+k-1-s-e
        if n+k-1-e <= s: return ''
        if s >= k-1: return self.seq[s-k+1:n-e]
        if e < n:  
            seq = self.seq[:n-e]
            e=n
        if e >= k-1: return compl(get_fasta(self.twin, k, e, s)[::-1]) + seq
    
        ### not whole sequence available!
        
        if s < n: seq_pocz = compl(get_fasta(self.twin, k, k-1, s)[::-1])
        else: seq_pocz = ''
        
        s = max(s, n) 
        
        #try without recursion
        for node in self.prev_nodes: 
            if node.length >= 2*(k-1) or node.n >= k-1-s  :
                missing = get_fasta_back(node, k, node.length-(k-1-s), e-n)
                return seq_pocz + missing + seq
        
        for node in self.next_nodes: 
            if node.length >= 2*(k-1) or node.n >= k-1-e :
                missing = get_fasta_forward(node, k, s-n, node.length-(k-1-e))
                return seq_pocz + missing + seq

        #with recursion
        
        missing = ''
        for node in self.prev_nodes: 
                newmissing = node.seq[-(k-1-s):node.n-(e-n)]
                if len(newmissing) > len(missing):
                    missing = newmissing
                    
        e += len(missing)
        seq = missing + seq
        
    
        for node in self.prev_nodes: 
            if node != self:
                print 'rec', node.id, node.length, node.length-(k-1-s), e-n
                missing =  get_fasta_back(node, k, node.length-(k-1-s), e-n)
                if missing != None: return seq_pocz + missing + seq
        return None
        




