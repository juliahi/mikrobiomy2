import SGAFilter
import cPickle
import time
import graph_stats
import sys


from common_parameters import *
from common import *



def run_stats(sg):
    count_sums1 = graph_stats.analyze_counts(sg)
    degrees1 = graph_stats.analyze_degrees(sg)
    sizes1 = graph_stats.analyze_connected_components(sg)
    lengths1 = graph_stats.analyze_lengths(sg)
    lfcs1 = graph_stats.analyze_foldchanges(sg)
    return [count_sums1, degrees1, sizes1, lengths1, lfcs1]


def load_graph(filename, folder, suf, conds):
    print "Loading graph:", filename, give_time()
    sg = SGAFilter.SgaGraph(conds)
    sg.init_graph(filename)
    print give_time()
    sg.load_graph()
    print "Finished loading graph", give_time()

    sg.add_duplicates_asqg(folder + "/merged" + suf + ".dups.merged" + suf + ".asqg")
    sg.add_duplicates_fasta(folder + "/merged" + suf + ".fa", reverse=True)
    sg.add_duplicates_fasta(folder + "/control" + suf + ".fa")
    sg.add_duplicates_fasta(folder + "/treated" + suf + ".fa")
    sg.finish_loading_counts()
    print "Finished loading duplicates", give_time()
    sys.stdout.flush()

    return sg



def simplify_graph(sg):
    stats = []
    sg.compress_simple_paths()
    sg.remove_short_islands(MIN_LENGTH)
    print "Finished simplifying paths", give_time()

    sys.stdout.flush()
    graph_stats.short_summary(sg)
    stats += run_stats(sg)
    print "---------", give_time(), "---------"
    sys.stdout.flush()

    # REMOVE DEAD-ENDS and simplify
    for i in xrange(DEADENDS_REMOVE_ROUNDS):
        sg.remove_deadends_by_length(DEADENDS_MIN_LENGTH)
        graph_stats.short_summary(sg)

        stats += run_stats(sg)
        print "Finished dead-ends round %d" % i, give_time()

    sg.compress_simple_paths()
    sg.remove_short_islands(MIN_LENGTH)

    print "Finished dead-ends removal", give_time()
    sys.stdout.flush()

    graph_stats.short_summary(sg)
    stats += run_stats(sg)

    print "---------", give_time(), "---------"
    sys.stdout.flush()
    return stats


def load_and_simplify():
    stats = []
    sg = load_graph(graph_filename, folder, suf, conditions)

    ### Write basic graph statistics:
    graph_stats.short_summary(sg)
    stats += run_stats(sg)
    print "---------", give_time(), "---------"
    sys.stdout.flush()

    ### simplify graph and compute graph statistics
    stats += simplify_graph(sg)

    ### save graph in ASQG and pickle, and stats as pickle
    sg.write_to_asqg(my_simplified_graph_filename, contigs=None, rename=False)
    cPickle.dump(stats, open(simplify_stats_filename, 'wb'), protocol=2)
    cPickle.dump(sg, open(my_renamed_graph_pickle, 'wb'), protocol=2)

    print "Finished saving stats", give_time()
    return sg


if __name__ == "__main__":

    #sg = load_and_simplify()

    # tmp: remove TODO
    sg = load_graph(my_simplified_graph_filename, folder, suf, conditions)
    sg.filename = my_simplified_graph_filename

    ### normalize graph and rename
    sg.normalize_counts()
    # writing to ASQG and renaming permanently:
    sg.write_to_asqg(my_renamed_graph_filename, contigs=my_renamed_graph_filename+"-contigs.fa",
                     rename=renamed_nodenames)

    ### save normalized and renamed graph in pickle
    cPickle.dump(sg, open(my_renamed_graph_pickle, 'wb'), protocol=2)
    print "Finished saving graph", give_time()


