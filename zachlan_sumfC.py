

import sys
from Graph import *


def find_important(vg): ##update?
    
        selectedUp = []
        selectedDown = []
        assembly_id = 1
        while True:
            max_node = vg.max_reads_node(assembly_id)
            if max_node == None:
                break
            if max_node.get_score() > 2:
                selectedUp.append(max_node.id)
            elif max_node.get_score() < 0.5:
                selectedDown.append(max_node.id)
            assembly_id += 1
            
        return selectedUp, selectedDown

"""    def traverse(vg, min_length): #travers wg najwyzszej liczby zliczen #TODO:update?
        selected = []
        assembly_id = 1
        while True:
            max_node = vg.max_reads_node(assembly_id)
            if max_node == None:
                break
            print max_node
            tmp_node = max_node
            current = []
            current.append(tmp_node)
            sum_score = node.get_score()
            while True:
                
                next_node = vg.get_max_next(tmp_node, assembly_id)
                print tmp_node, tmp_node.get_score()
                if next_node == None:
                    break
                tmp_node = next_node
                current.append(tmp_node)
                sum_score += tmp_node.get_score()
            
            tmp_node = max_node
            while True:
                
                next_node = vg.get_max_prev(tmp_node, assembly_id)
                print tmp_node, tmp_node.get_score()
                if next_node == None:
                    break
                tmp_node = next_node 
                current = [tmp_node] + current
                sum_score += tmp_node.get_score()
            
            #"usuwamy" srednie pokrycie
            
            print "select:", assembly_id
            for node in current:
                print node
                node.select(sum_score/len(current))
            selected.append((current, sum_score))
            assembly_id += 1
        return selected
"""


def find_diff(vg, min_length, min_fc):
    assemblies = []
    assembly_id = 1
    vg.start_queue()
    while True:
        node = vg.max_fc_node(assembly_id)
        if not node: break
        path = Path(node, assembly_id)
        while path.expand(min_fc):
            pass
        fc = path.foldChange()
        if path.length >= min_length and (fc >= min_fc or 1/fc >= min_fc):
            path.select()
            assemblies.append(path)
        assembly_id += 1
    return assemblies


def write_paths(vg, paths, output_name):
    with open(output_name+"assemblies.txt", 'w') as out_txt:
        with open(output_name+"assemblies.fa", 'w') as out_fa:
            out_txt.write("ID\tfoldChange\tlength\tnodes\tnodeFC\n")
            for i, path in enumerate(paths):
                #if len(path.node_ids()) < 2: continue
                nodes = path.get_nodes()
                fa = vg.get_fasta(nodes)
                fcs = map(lambda node: str(node.foldChange()), nodes)
                
                #print path.length, len(fa)
                out_txt.write("%d\t%f\t%f\t%s\t%s\n"%(path.assembly_id, path.foldChange(), len(fa), 
                                                  ','.join(map(str, path.node_ids())),
                                                  ','.join(fcs)))
                #out_fa.write('>%s\n%s\n%s\n%s\n'%(path.name(), fa, compl(fa), compl(fa[::-1])))
                out_fa.write('>%s\n%s\n'%(path.name(), fa))
                
        

def count_reads(fq_file):
    #def blocks(files, size=65536):
        #while True:
            #b = files.read(size)
            #if not b: break
            #yield b

    #with open(fq_file, "r") as f:
        #return sum(bl.count("\n") for bl in blocks(f))/4
    suma=0
    for line in open(fq_file, 'r'):
        suma += 1
    return suma/4


import argparse
if __name__ == "__main__":
    
    def pair(arg):
        return map(int, arg.split(','))
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--files', type=pair, nargs='+', help='for each used sample (fq.1.gz) give a pair isPositive,n.of_reads')
    parser.add_argument('-g', '--graph', help='LastGraph')
    parser.add_argument('-o', '--output', help='output name')
    parser.add_argument('--paired', action='store_true', help='where fq files paired')
    
    #if len(sys.argv) <= 4 or len(sys.argv[3].split(',')) != len(sys.argv[4].split(',')):
        #print """Usage:
        #Graph.py LastGraph -output_name file1.fq,file2.fq[,file3.fq...] cond1,cond2[,cond3...]
        #where condN is 0 or 1
        #"""
    
    args = parser.parse_args()
   
    print args.files
    reads_in_files = [x[1] for x in args.files]
    conds = [x[0] for x in args.files]
    
    if args.paired:
        reads_in_files = [2*x for x in reads_in_files]
    
    vg = VelvetGraph(args.graph, reads_in_files, conds)
    
    #MINLEN = 2*vg.k                #TODO: as program parameters?
    MINLEN = 200                
    MINFC = 2
    
    
    #for i in xrange(1,vg.n_nodes+1):
        ##f = vg.get_fasta_ids([i])
        #f2 = vg.get_fasta_ids([-i])
        
        #print i, #f==compl(f2)[::-1]
        ##print f
        #print f2
    
    
    ####up, down = find_important(vg)
    ####assemblies = traverse(vg)
    
    
    paths = find_diff(vg, MINLEN, MINFC)
    
    write_paths(vg, paths, args.output)
    
