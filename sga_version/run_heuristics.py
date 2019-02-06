#import SGAFilter
import cPickle
from load_filter_graph import give_time
import heuristics
from common import *

# test_name = "sga_test_notrim"
# conds = {'6685_04-06-2015': 0, '6685_16-06-2015': 1}
# folder = "/mnt/chr7/data/julia/"+test_name
# suf = ".preprocessed_q20.ec.filter.pass.rmdup"
# filename = folder+"/merged"+suf+"_51.asqg"


conds = {'6683_16-06-2015':1, '6685_04-06-2015':0, '6685_16-06-2015':1,
         '6690_04-06-2015':0, '6690_16-06-2015':1, '6695_04-06-2015':0,
         '6695_16-06-2015':1, '6704_04-06-2015':0, '6704_16-06-2015':1}
test_name = 'sga_test_full_notrim_paired_reversed'
folder = "/mnt/chr7/data/julia/"+test_name
suf = ".preprocessed_qf5.ec.filter.pass.rmdup"
filename = folder+"/merged"+suf+"_31.asqg"


MIN_LENGTH = 500
MIN_FC = 2

print "Loading starts", give_time()

# sg = cPickle.load(open('simplified_'+test_name+suf+'_sga.pickle', 'wb'), protocol=2)
# print "Finished loading graph", give_time()
# sg1 = SGAFilter.get_random_components(sg, 1)

# cPickle.dump(comp, open('1randomcomponent_'+test_name+suf+'_sga.pickle', 'wb'), protocol=2)
# print "Finished dumping component", give_time()
sg1 = cPickle.load(open('1randomcomponent_'+test_name+suf+'_sga.pickle', 'rb'))
print "Finished loading graph", give_time()


outname = "simplified_1random"
# Heuristic 1: take longest
print "Start heuristic 1: longest", give_time()
longest = heuristics.take_longest(sg1)
print "Found %d paths" % len(longest)
# filter by FC
longest_filter = sg1.save_path_sequences(longest, outname+'_sg1_longest.fa', min_fc = MIN_FC)
print "Found %d paths with foldchange >= %f" % (len(longest_filter), MIN_FC)
print "Finished heuristic 1: longest", give_time()

# Heuristic 2: take longest until foldchange
print "Start heuristic 2: longest until FC", give_time()
longest_fc = heuristics.take_longest_minfc(sg1, MIN_FC)
print "Found %d paths with foldchange >= %f" % (len(longest_fc), MIN_FC)
sg1.save_path_sequences(longest_fc, outname+'_sg1_longestfc.fa')
print "Finished heuristic 2: longest until FC", give_time()

# Heuristic 3: take best foldchange until > FC
print "Start heuristic 3: take best foldchange until > FC", give_time()
best_fc = heuristics.take_best_fc(sg1, MIN_FC)
print "Found %d paths with foldchange >= %f" % (len(best_fc), MIN_FC)
sg1.save_path_sequences(best_fc, outname+'_sg1_bestfc.fa')
print "Finished heuristic 3: take best foldchange until > FC", give_time()









