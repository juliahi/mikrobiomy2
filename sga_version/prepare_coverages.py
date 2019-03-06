from common import *
from common_parameters import *
from compare_assemblies import *
import glob
import sys
import cPickle


if __name__ == "__main__":

    if sys.argv[1] == "singleend":
        workdir = OUTPUTDIR+'heuristics_' + outname + '/'
        mapping_sufix = "/singleend_newkallisto/*_sorted.bam"
    elif sys.argv[1] == "bowtie":
        workdir = OUTPUTDIR+'heuristics_' + outname + '/'
        mapping_sufix = "/bowtie_results/*_sorted.bam"
    else:
        sys.exit(1)

    names = []
    seqs = []
    counts = []
    sums = {}
    for name in ["longestfc"]:
        bam_files = glob.glob(workdir + name + mapping_sufix)
        tmp_names, tmp_seqs, tmp_counts, tmp_sums = init_assembly(workdir + name + "_filtered_fc2.fa", bam_files)
        names.append(tmp_names)
        seqs.append(tmp_seqs)
        counts.append(tmp_counts)
        sums[name] = tmp_sums

    sga_names, sga_seqs, sga_counts, sga_sums = init_assembly(OUTPUTDIR + 'sga_scaffold_' + outname + '/sga-scaffolds_200.fa',
                                                              glob.glob(OUTPUTDIR + 'sga_scaffold_' + outname + mapping_sufix))

    m_names, m_seqs, m_counts, m_sums = init_assembly(OUTPUTDIR+'megahit/megahit_200.fa',
                                                      glob.glob(OUTPUTDIR+"megahit"+mapping_sufix))
    names += [sga_names, m_names]
    seqs += [sga_seqs, m_seqs]
    counts += [sga_counts, m_counts]
    sums["sga"] = sga_sums
    sums["megahit"] = m_sums

    cPickle.dump([names, seqs, counts, sums], open(sys.argv[1]+'_filtered.pickle', 'wb'), protocol=2)
