
from read_kallisto_stats import *

import sys, glob

#dirname="/mnt/chr7/data/julia/kallisto_stats/"
dirname=sys.argv[1]
outname=sys.argv[2]
kvals=sys.argv[3:]


def mean(l): return sum(l)/len(l)





with open(outname+'.tsv', 'w') as outfile:
    with open(outname+'.tex', 'w') as outfile2:
        outfile.write("Filename\tmapped\tuniquelly mapped\tno kmers\tconflicts\tmean number of kmers\n")
        outfile2.write("\\begin{tabular}{l|rrrr} \n")
        outfile2.write("Filename & mapped (uniquelly) & no kmers & conflicts & mean number of kmers \\\\ \n")
        for k in kvals:
            k = int(k) 
            
            mapowalne=[]
            umapowalne=[]
            niemapowalne_brak=[]
            niemapowalne_konflikt=[]
            kmers=[]

            outfile2.write("\\hline \\\\\n")

            #name="6685_04-06-2015_depl"
            print(sorted(glob.glob(dirname+("/*_kallisto_%d_out/stats.txt"%k))))

            for filename in sorted(glob.glob(dirname+("/*_kallisto_%d_out/stats.txt"%k))):
                
                stats=Stats(filename)
                print filename, sum(stats.mapped()), sum(stats.unique())

                mapowalne.append(stats.fraction_mapped())
                umapowalne.append(stats.fraction_unique())
                kmers.append(sum(stats.data["kmers1"] + stats.data["kmers2"])/stats.n)
                niemapowalne_brak.append(stats.fraction_zerokmers())
                niemapowalne_konflikt.append(stats.fraction_conflicts())
                
                
                outfile.write("%s\t%f\t%f\t%f\t%f\t%f\n"%(filename, mapowalne[-1], umapowalne[-1], niemapowalne_brak[-1], niemapowalne_konflikt[-1], kmers[-1]))
                outfile2.write("%s, k=%d & %.1f \\%% (%.1f \\%%) & %.1f \\%% & %.1f \\%% & %.3f  \\\\ \n"%(filename.split('/')[-2][:15].replace('_', '\_'), k, 
                    mapowalne[-1]*100, umapowalne[-1]*100, niemapowalne_brak[-1]*100, niemapowalne_konflikt[-1]*100, kmers[-1]))
            if mapowalne != []:
                outfile.write("kmer=%d mean\t%f\t%f\t%f\t%f\t%f\n"%(k, mean(mapowalne), mean(umapowalne), mean(niemapowalne_brak), mean(niemapowalne_konflikt), mean(kmers)))
                outfile2.write("k=%d mean & %.1f \\%% (%.1f \\%%) & %.1f \\%% & %.1f \\%% & %.3f \\\\ \n"%(k, 
                    mean(mapowalne)*100, mean(umapowalne)*100, mean(niemapowalne_brak)*100, mean(niemapowalne_konflikt)*100, mean(kmers)))
            
        outfile2.write("\\end{tabular} \n")


