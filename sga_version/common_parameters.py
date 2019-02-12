### parameters of experiment and file paths

conditions = {'6683_16-06-2015': 1, '6685_04-06-2015': 0, '6685_16-06-2015': 1,
             '6690_04-06-2015': 0, '6690_16-06-2015': 1, '6695_04-06-2015': 0,
             '6695_16-06-2015': 1, '6704_04-06-2015': 0, '6704_16-06-2015': 1}
test_name = 'sga_test_full_notrim_paired_reversed'
folder = "/mnt/chr7/data/julia/" + test_name
suf = ".preprocessed_qf5.ec.filter.pass.rmdup"

K=31

graph_filename = folder + "/merged" + suf + "_" + str(K) + ".asqg"

# simplification and filtering
MIN_LENGTH = 200
DEADENDS_REMOVE_ROUNDS = 1
DEADENDS_MIN_LENGTH = 200

OUTPUTDIR="/home/julia/mikrobiomy_results/"
#my_simplified_graph_filename = OUTPUTDIR + "my_simplified" + str(DEADENDS_REMOVE_ROUNDS) + "_merged" + suf + "_" + str(K) +  "_"+str(MIN_LENGTH)+".asqg"
outname = "mysimplified%d_%d"%(DEADENDS_REMOVE_ROUNDS, MIN_LENGTH)
my_simplified_graph_filename = "%s%s.asqg" % (OUTPUTDIR, outname)
my_simplified_graph_pickle = "%s%s.pickle" % (OUTPUTDIR, outname)
renamed_nodenames = "%s%s_nodes_renamed.tsv" % (OUTPUTDIR, outname)
my_renamed_graph_pickle = "%s%s_renamed.pickle" % (OUTPUTDIR, outname)
my_renamed_graph_filename = "%s%s_renamed.asqg" % (OUTPUTDIR, outname)

simplify_stats_filename = "%s%s_stats.pickle" % (OUTPUTDIR, outname)



# for plots
import seaborn as sns
color_palette = sns.color_palette('muted6')
labels = ["max_len", "max_len_fc", "best_fc", "sga", "megahit"]
colors = {"max_len": color_palette[0], "max_len_fc": color_palette[5], "best_fc": color_palette[2],  \
          "sga": color_palette[4], "megahit": color_palette[1]}

labels = ["longest", "longestfc", "bestfc", "sga", "megahit"]
colors = {"longest": color_palette[0], "longestfc": color_palette[5], "bestfc": color_palette[2],  \
          "sga": color_palette[4], "megahit": color_palette[1]}


# for heuristics

MIN_FC = 2

longest_fasta = "%slongest_filtered_fc%d.fa" % (OUTPUTDIR, MIN_FC)
longestfc_fasta = "%slongestfc_filtered_fc%d.fa" % (OUTPUTDIR, MIN_FC)
bestfc_fasta = "%sbestfc_filtered_fc%d.fa" % (OUTPUTDIR, MIN_FC)
sgafile = "%ssga-scaffolds.fa" % OUTPUTDIR
sga_fasta = "%ssga_filtered_fc%d.fa" % (OUTPUTDIR, MIN_FC)


NPAIRS = 90582474


