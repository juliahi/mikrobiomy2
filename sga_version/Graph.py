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
    epsilon=.001
    def __init__(self, id, seq):
        self.id = id
        self.next_edges = [] 
        self.prev_edges = [] 
        self.nreads = [Node.epsilon, Node.epsilon] #defaultdict(0.)
        #self.n = int(values[0]) ##=len(seq)
        
        self.selected = 0  #assembly number in which it was selected lastly
        ####self.seq = seq
        self.length = len(seq)
        

    def foldChange(self):
        if self.nreads[1] == 0: return None
        return 1.*self.nreads[0] / self.nreads[1]
    
    
    def get_score(self):    # should return float!
        fc = self.foldChange()
        if fc > 1: fc = 1./fc
        return fc
        
    def add_edge(self, edge, first):
        if first == 0:
            self.next_edges.append(edge)
        else:
            self.prev_edges.append(edge)
    
    def neighbours(self): #TODO: co jesli zmienia sie kierunek?
        return [e.node2 for e in self.next_edges] + [e.node1 for e in self.prev_edges]
    


class Edge:
    def __init__(self, values, node1, node2):
        #self.name1 = values[0]
        self.start1 = int(values[2])
        self.end1 = int(values[3])
        #self.len1 = int(values[4])
        #self.name2 = values[1]
        self.start2 = int(values[5])
        self.end2 = int(values[6])
        #self.len2 = int(values[7])
        
        node1.add_edge(self, 0)
        node2.add_edge(self, 1)
        self.node1 = node1
        self.node2 = node2
        
        self.orient = int(values[8]) #(0=same, 1=reverse) 
        #self.diffs = int(values[9]) 
        
    def whattype(self):
        if self.node1.length == abs(self.end1 - self.start1)+1: return 'incl1'
        if self.node2.length == abs(self.end2 - self.start2)+1: return 'incl2'
        if self.orient == 0: 
            if self.start2 < self.end2: return '>->'
            else: return '<-<'
        else:
            if self.start2 < self.end2: return '>-<'
            else: return '<->'

    def join_nodes(self):
        if len(self.node1.next_edges) == 1 and len(self.node2.prev_edges) == 1: #join
            newid = self.node1.id + '|' + self.node2.id
            if self.orient == 0:
                #newseq = self.node1.seq += self.node2.seq[self.end2+1:]
                self.node1.next_edges = self.node2.next_edges
                x = self.start1
                for e in self.node1.next_edges:
                    e.start1 += x
                    e.end1 += x
                    e.node1 = self.node1
                self.node1.id = newid
                self.node1.nreads = [c[0] + c[1] for c in zip(self.node1.nreads, self.node2.nreads)]
                self.node1.length += self.end2 - self.start2 + 1
                return True
        return False
            

class SgaGraph:
    
    def __init__(self, filename, conds=None): #1 or 2 conditions
        self.conds = conds
        self.nodes = {}
        self.edges = []
        self.duplicates_dict = {}
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
                while line: #Edges
                    self.add_edge(line.split()[1:])
                    line = f.readline()
            except Exception as e: 
                print e
                
                pass
            print 'no. edges loaded:', len(self.edges)
    
    
    def whatcond(self, node_id):
        cond = node_id.split(':')[-1]
        if cond[-2] == '/': cond=cond[:-2]
        #TODO: kontrola bledu?
        return self.conds[cond]
    
    def node(self, node_id):
        try: return self.nodes[node_id]
        except: 
            try: return self.nodes[self.duplicates_dict[node_id]]
            except: return None
            

    def add_node(self, node_id, seq):
        node = Node(node_id, seq)
        if self.conds:
            node.nreads[self.whatcond(node_id)] += 1
        self.nodes[node_id] = node
        
    def add_edge(self, values):
        e = Edge(values, self.node(values[0]), self.node(values[1]))
        self.edges.append(e)  #?
        
    def add_dupl(self, name, count, reverse):
        node = self.node(name)
        if node:
            if reverse:
                node.nreads[1-self.whatcond(name)] += count
            else:
                node.nreads[self.whatcond(name)] += count
            
                                    
    def add_duplicates_fasta(self, filename, reverse=False): 
        #reverse - add to the second condition because duplicate is from another condition
        #filename = fasta file after rmdup (not duplicated) -- contains counts 
        with open(filename) as f:
            try:
                added=0
                addedc=0
                while True: 
                    line = f.readline()
                    if line and line[0] == '>':
                        _, name, number = line.split()
                        c = int(number.split('=')[1])-1
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
            print 'added %d reads to %d sequences'%(addedc, added)
        
    def add_duplicates_asqg(self, filename): 
        # filename = asqg file of overlap btw duplicates 
        # create dictionary of removed reads: existing reads
        with open(filename) as f:
            f.readline()
            try:
                while True: #Nodes
                    line = f.readline()
                    if line[0] != 'V':
                        print line
                        break
                #lines = 0
                countdups = 0
                while True: #Edges
                    if line and line[0] == 'E':
                        line = line.split()[1:]
                        #lines += 1
                        #if lines % 100000 == 0: print lines
                        name1 = line[0].split(',')[0]           # remove seqrank=X from name of duplicated seq
                        if (line[2] == '0' and int(line[3]) == int(line[4]) - 1) or (line[5] == '0' and int(line[6]) == int(line[7]) - 1):
                            if self.whatcond(name1) != self.whatcond(line[1]): 
                                countdups += 1
                                self.duplicates_dict[name1] = line[1] 
                                
                    else: break    
                    line = f.readline()
            except Exception as e: 
                print e, 'line:', line, "."
                pass
            print "Found %d duplicates"%countdups
    
    def compress_simple_paths(self):
        print "Edges before compression:", len(self.edges)
        
        newedges = []
        for edge in self.edges:
            oldname1 = edge.node1.id
            oldname2 = edge.node2.id
            
            if not edge.join_nodes(): #not joined
                newedges.append(edge)
            else: #joined, remove node
                self.nodes[edge.node1.id] = self.nodes.pop(oldname1)
                del self.nodes[oldname2]
        self.edges = newedges
        print "After compression:", len(self.edges)
    
    def node_outdegrees(self):
        return [len(x.next_edges)  for x in self.nodes.values()]
    def node_indegrees(self):
        return [len(x.prev_edges)  for x in self.nodes.values()]
    def edge_types(self):
        return [e.whattype()  for e in self.edges]
    
    def traverse(self, nodename, visited):
        q = Q.Queue() # queue with nodes
        #if not visited[nodename]:
        q.put(self.nodes[nodename])
        visited[nodename] = True
        size = 1
        while not q.empty():
            node = q.get()
            for neighbour in node.neighbours():
                if not visited[neighbour.id]:
                    visited[neighbour.id] = True
                    q.put(neighbour)
                    size += 1
                    
        return size
        
    
    def connected_components(self):
        visited = {x:False for x in self.nodes.keys()}
        components = 0
        sizes = []
        for name in visited.keys():
            if visited[name]:
                continue
            else:
                components += 1
                size = self.traverse(name, visited)
                sizes.append(size)
                if components > 1000000: break
        return components, sizes
    
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
    
    def remove_edge(self, edge):
        try: edge.node1.next_edges.remove(edge)
        except: edge.node1.prev_edges.remove(edge)
        
        try: edge.node2.prev_edges.remove(edge)
        except: edge.node2.next_edges.remove(edge)
        
    
    def clean_edges(self):
        ##removes edges with >-< orientation
        for edge in self.edges:
            if edge.orient == 1: self.remove_edge(edge)
        self.edges = [e for e in self.edges if e.orient == 0]    





############do testowania
#def get_fasta(selfnode, k, s=0, e=0): #s,e - offsets from start/end, 0-based
        #seq = ''
        #self = selfnode
        #n = self.n
        #print self.id,s,e, n+k-1-s-e
        #if n+k-1-e <= s: return ''
        #if s >= k-1: return self.seq[s-k+1:n-e]
        #if e < n:  
            #seq = self.seq[:n-e]
            #e=n
        #if e >= k-1: return compl(get_fasta(self.twin, k, e, s)[::-1]) + seq
    
        #### not whole sequence available!
        
        #if s < n: seq_pocz = compl(get_fasta(self.twin, k, k-1, s)[::-1])
        #else: seq_pocz = ''
        
        #s = max(s, n) 
        
        ##try without recursion
        #for node in self.prev_edges: 
            #if node.length >= 2*(k-1) or node.n >= k-1-s  :
                #missing = get_fasta(node, k, node.length-(k-1-s), e-n)
                #return seq_pocz + missing + seq
        
        #for node in self.next_edges: 
            #if node.length >= 2*(k-1) or node.n >= k-1-e :
                #missing = get_fasta(node, k, s-n, node.length-(k-1-e))
                #return seq_pocz + missing + seq

        ##with recursion
        
        #print self.id, seq_pocz, seq, s, e
        
        #missing = ''
        #for node in self.prev_edges: 
                #newmissing = node.seq[-(k-1-s):node.n-(e-n)]
                #if len(newmissing) > len(missing):
                    #missing = newmissing
                    
        #e += len(missing)
        #seq = missing + seq
        #missing = ''
        #for node in self.next_edges:
                #newmissing = get_fasta(node, k, s-n, k-1)  #only from twin seq
                #if newmissing != None and len(newmissing) > len(missing):
                    #missing = newmissing
        #s += len(missing)
        #seq_pocz += missing
        
        #print self.id, seq_pocz, seq
    
        #for node in self.prev_edges: 
            #if node != self:
                #print 'rec', node.id, node.length, node.length-(k-1-s), e-n
                ##missing =  'N'*(node.length-e+n-node.length+(k-1-s)) #get_fasta(node, (k, node.length-(k-1-s), e-n)
                #missing =  get_fasta_back(node, k, node.length-(k-1-s), e-n)
                #if missing != None: return seq_pocz + missing + seq
        #for node in self.next_edges: 
            #if node != self:
                #print 'rec2', node.id, node.length, s-n, node.length-(k-1-e)
                ##missing =  'N'*(node.length-s+n-node.length+(k-1-e))  ###get_fasta(node, (k, s-n, node.length-(k-1-e))
                #missing =  get_fasta_forward(node, k, s-n, node.length-(k-1-e))
                #if missing != None: return seq_pocz + missing + seq
        #return None
        




#def get_fasta_forward(selfnode, k, s=0, e=0): #s,e - offsets from start/end, 0-based
        #seq = ''
        #self=selfnode
        #n = self.n
        #if n+k-1-e <= s: return ''
        #if s >= k-1: return self.seq[s-k+1:n-e]
        #if e < n:  
            #seq = self.seq[:n-e]
            #e=n
        #if e >= k-1: return compl(get_fasta(self.twin, k, e, s)[::-1]) + seq
    
        #### not whole sequence available!
        
        #if s < n: seq_pocz = compl(get_fasta(self.twin, k, k-1, s)[::-1])
        #else: seq_pocz = ''
        
        #s = max(s, n) 
        #for node in self.prev_edges: 
            #if node.length >= 2*(k-1) or node.n >= k-1-s  :
                #missing = get_fasta_back(node, k, node.length-(k-1-s), e-n)
                #return seq_pocz + missing + seq
        #for node in self.next_edges: 
            #if node.length >= 2*(k-1) or node.n >= k-1-e :
                #missing = get_fasta_forward(node, k, s-n, node.length-(k-1-e))
                #return seq_pocz + missing + seq

        ##with recursion
        #missing = ''
        #for node in self.next_edges:
                #newmissing = get_fasta_forward(node, k, s-n, k-1)  #only from twin seq
                #if newmissing != None and len(newmissing) > len(missing):
                    #missing = newmissing
        #s += len(missing)
        #seq_pocz += missing
        
        #for node in self.next_edges: 
            #if node != self:
                #missing = get_fasta_forward(node, k, s-n, node.length-(k-1-e))
                #if missing != None: return seq_pocz + missing + seq
        #return None
        
        

#def get_fasta_back(selfnode, k, s=0, e=0): #s,e - offsets from start/end, 0-based
        #seq = ''
        #self = selfnode
        #n = self.n
        #print self.id,s,e, n+k-1-s-e
        #if n+k-1-e <= s: return ''
        #if s >= k-1: return self.seq[s-k+1:n-e]
        #if e < n:  
            #seq = self.seq[:n-e]
            #e=n
        #if e >= k-1: return compl(get_fasta(self.twin, k, e, s)[::-1]) + seq
    
        #### not whole sequence available!
        
        #if s < n: seq_pocz = compl(get_fasta(self.twin, k, k-1, s)[::-1])
        #else: seq_pocz = ''
        
        #s = max(s, n) 
        
        ##try without recursion
        #for node in self.prev_edges: 
            #if node.length >= 2*(k-1) or node.n >= k-1-s  :
                #missing = get_fasta_back(node, k, node.length-(k-1-s), e-n)
                #return seq_pocz + missing + seq
        
        #for node in self.next_edges: 
            #if node.length >= 2*(k-1) or node.n >= k-1-e :
                #missing = get_fasta_forward(node, k, s-n, node.length-(k-1-e))
                #return seq_pocz + missing + seq

        ##with recursion
        
        #missing = ''
        #for node in self.prev_edges: 
                #newmissing = node.seq[-(k-1-s):node.n-(e-n)]
                #if len(newmissing) > len(missing):
                    #missing = newmissing
                    
        #e += len(missing)
        #seq = missing + seq
        
    
        #for node in self.prev_edges: 
            #if node != self:
                #print 'rec', node.id, node.length, node.length-(k-1-s), e-n
                #missing =  get_fasta_back(node, k, node.length-(k-1-s), e-n)
                #if missing != None: return seq_pocz + missing + seq
        #return None
        




