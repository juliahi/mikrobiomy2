
OUTDIR=/mnt/chr7/data/julia/sga_test

sh run_sga_separate.sh
exit
python graph_stats.py -g $OUTDIR/*51.asqg -o $OUTDIR/out

python graph_stats.py -g $OUTDIR/*51-graph.asqg -o $OUTDIR/out_simplified
