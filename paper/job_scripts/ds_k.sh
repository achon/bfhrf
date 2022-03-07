n=100
for k in 1000 25000 50000 75000 100000
do
  echo $k : $(date)
  command time -f "CPU:%P  MEM(kbytes):%M" python3 DendropySingle.py gene_trees.1.99.$k gene_trees.1.99.$k $n &> gene_trees.1.99.$k.ds.log
done
