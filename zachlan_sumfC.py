

import sys
from Graph import *
import matplotlib
matplotlib.use('Agg')
from datetime import datetime

import math


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
    suma=0
    for line in open(fq_file, 'r'):
        suma += 1
    return suma/4


import seaborn as sns
import pandas

def plot_fc(vg, outputdir):
    df = pandas.DataFrame.from_dict({'mean_counts': [np.mean(node.nreads) for node in vg.nodes[::2]], 
                                     'log2foldChange':[math.log(node.foldChange(),2) for node in vg.nodes[::2]]  })
    plot = sns.jointplot(x="mean_counts", y="log2foldChange", data=df)
    plot.savefig(outputdir+'_fc_nodes.pdf')
        



import argparse    
from cProfile import Profile
from pstats import Stats
prof = Profile()
prof.disable() 


if __name__ == "__main__":
    def pair(arg):
        return map(int, arg.split(','))
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--files', type=pair, nargs='+', help='for each used sample (fq.1.gz) give a pair isPositive,n.of_reads')
    parser.add_argument('-g', '--graph', help='LastGraph')
    parser.add_argument('-o', '--output', help='output name')
    parser.add_argument('--paired', action='store_true', help='where fq files paired')
    parser.add_argument('--profile', action='store_true', help='cProfile')
    parser.add_argument('--minlen', type=int, default=None, help='minimum length of transcript, default=2*k')
    parser.add_argument('--minfc', type=float, default=2., help='minimum foldChange of transcript, default=2')
    
    args = parser.parse_args()
   
    reads_in_files = [x[1] for x in args.files]
    conds = [x[0] for x in args.files]
    
    if args.paired:
        reads_in_files = [2*x for x in reads_in_files]
    
    
    if args.profile:
        prof.enable()
    
    with open(args.output+'.log', 'w+') as log:
        log.write("Running with parameters: %s\n\n"%(str(args)))
        
        vg = VelvetGraph(args.graph, reads_in_files, conds)
    
        if args.minlen == None:
            args.minlen = 2*vg.k
        log.write("%s: Graph loaded. Minlen: %d\n"%(datetime.now(), args.minlen))
        
        
        if args.profile:
            prof.disable()
        plot_fc(vg, args.output)
        
        if args.profile:
            prof.enable()
        
        log.write("%s: Finding paths...\n"%datetime.now(), )
        log.flush()
        paths = find_diff(vg, args.minlen, args.minfc)
        
        log.write("%s: Found %d paths\n"%(datetime.now(), len(paths)))
        log.flush()
        
        write_paths(vg, paths, args.output)
    
    if args.profile:
        prof.disable()  # don't profile the generation of stats
        prof.dump_stats(args.output+'profile.stats')
        with open(args.output+'profile_stats.txt', 'a+') as output:
            stats = Stats(args.output+'profile.stats', stream=output)
            stats.sort_stats('cumulative', 'time')
            stats.print_stats()

    
    ####up, down = find_important(vg)
    ####assemblies = traverse(vg)
    
    
    
