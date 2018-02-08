
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



import Graph

import argparse    
import seaborn 


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--graph', help='ASQG file')
    parser.add_argument('-o', '--output', help='output name')
    parser.add_argument('--profile', action='store_true', help='cProfile')
    #parser.add_argument('--minlen', type=int, default=None, help='minimum length of transcript, default=2*k')
    #parser.add_argument('--minfc', type=float, default=2., help='minimum foldChange of transcript, default=2')
    
    args = parser.parse_args()
    
    if args.profile:
        from cProfile import Profile
        from pstats import Stats
        prof = Profile()
        #prof.disable() 
        prof.enable()
    
    if not args.output:
        args.output=args.graph + '_output'
    
    sg = SgaGraph(args.graph)
    
    node_out = sg.node_outdegrees()
    hist = sns.distplot(node_out)
    hist.savefig(args.output+'.png')
    
    #with open(args.output+'.log', 'a+') as log:
        #log.write("%s: Running with parameters: %s\n\n"%(datetime.now(), str(args)))
        #log.flush()
                
        #vg = VelvetGraph(args.graph, reads_in_files, conds)
    
        #if args.minlen == None:
            #args.minlen = 2*vg.k
        #log.write("%s: Graph loaded. Minlen: %d\n"%(datetime.now(), args.minlen))
        #log.flush()
        
        #if args.profile: ##TMP
            #prof.disable()
        
        
        ##if args.profile:
            ##prof.disable()
        ##plot_fc(vg, args.output)
        
        ##if args.profile:
            ##prof.enable()
        

        #log.write("%s: Finding paths...\n"%datetime.now(), )
        #log.flush()
        
        
        #paths_count = find_diff(vg, args.minlen, args.minfc, args.output)
        
        #log.write("%s: Found %d paths\n"%(datetime.now(), paths_count))
        #log.flush()
        
    
    if args.profile:
        prof.disable()  # don't profile the generation of stats
        prof.dump_stats(args.output+'profile.stats')
        with open(args.output+'profile_stats.txt', 'a+') as output:
            stats = Stats(args.output+'profile.stats', stream=output)
            stats.sort_stats('cumulative', 'time')
            stats.print_stats()

    
    ####up, down = find_important(vg)
    ####assemblies = traverse(vg)
    

    
