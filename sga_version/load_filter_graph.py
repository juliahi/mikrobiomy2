import SGAFilter
import cPickle
import time
import graph_stats
import sys
#
# test_name = "sga_test_notrim"
# conds = {'6685_04-06-2015': 0, '6685_16-06-2015': 1}
# folder = "/mnt/chr7/data/julia/"+test_name
# suf = ".preprocessed_q20.ec.filter.pass.rmdup"
# filename = folder+"/merged"+suf+"_51.asqg"


def give_time():
    return time.asctime(time.localtime(time.time()))


def run_stats(sg):
    count_sums1 = graph_stats.analyze_counts(sg)
    degrees1 = graph_stats.analyze_degrees(sg)
    sizes1 = graph_stats.analyze_connected_components(sg)
    lengths1 = graph_stats.analyze_lengths(sg)
    lfcs1 = graph_stats.analyze_foldchanges(sg)
    return [count_sums1, degrees1, sizes1, lengths1, lfcs1]


if __name__ == "__main__":
    conds = {'6683_16-06-2015': 1, '6685_04-06-2015': 0, '6685_16-06-2015': 1,
             '6690_04-06-2015': 0, '6690_16-06-2015': 1, '6695_04-06-2015': 0,
             '6695_16-06-2015': 1, '6704_04-06-2015': 0, '6704_16-06-2015': 1}
    test_name = 'sga_test_full_notrim_paired_reversed'
    folder = "/mnt/chr7/data/julia/" + test_name
    suf = ".preprocessed_qf5.ec.filter.pass.rmdup"
    filename = folder + "/merged" + suf + "_31.asqg"

    MIN_LENGTH = 300

    stats = []

    print "Loading graph:", filename,  give_time()
    sg = SGAFilter.SgaGraph(conds)
    sg.init_graph(filename)
    print give_time()
    sg.load_graph()
    print "Finished loading graph", give_time()

    sg.add_duplicates_asqg(folder+"/merged"+suf+".dups.merged"+suf+".asqg")
    sg.counts = None
    sg.add_duplicates_fasta_todict(folder+"/merged"+suf+".fa", reverse=True)
    sg.add_duplicates_fasta_todict(folder+"/control"+suf+".fa")
    sg.add_duplicates_fasta_todict(folder+"/treated"+suf+".fa")
    print "Finished loading duplicates", give_time()
    sys.stdout.flush()

    # Write basic graph statistics:
    graph_stats.short_summary(sg)

    stats += run_stats(sg)

    print "---------", give_time(), "---------"
    sys.stdout.flush()


    # COLLAPSING SIMPLE PATHS
    sg.compress_simple_paths()
    sg.remove_short_islands(MIN_LENGTH)

    print "Finished simplification of paths", give_time()
    sys.stdout.flush()
    graph_stats.short_summary(sg)

    stats += run_stats(sg)

    print "---------", give_time(), "---------"
    sys.stdout.flush()

    # REMOVE DEAD-ENDS and simplify
    sg.remove_deadends_by_length(MIN_LENGTH)
    sg.compress_simple_paths()
    sg.remove_short_islands(MIN_LENGTH)

    print "Finished dead-ends removal", give_time()
    sys.stdout.flush()

    graph_stats.short_summary(sg)

    stats += run_stats(sg)

    print "---------", give_time(), "---------"
    sys.stdout.flush()


    # Removing simple nodes, saving graph
    # node_list, seqs = sg.save_simple_nodes(MIN_LENGTH, test_name+'_simpleseqs'+suf+'_sga.fa', filename)
    # print "Finished saving not connected nodes with length >= %d: %s" % (MIN_LENGTH, give_time())
    # graph_stats.short_summary(sg)

    stats += run_stats(sg)

    cPickle.dump(stats, open('stats_'+test_name+suf+'_sga.pickle', 'wb'), protocol=2)
    print "Finished saving stats", give_time()

    cPickle.dump(sg, open('simplified_'+test_name+suf+'_sga.pickle', 'wb'), protocol=2)

    print "Finished saving graph", give_time()


