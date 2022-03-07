k=1000
for n in 99 249 499 749 999
do
  ((m=n+1))
  echo "$m" : $(date)
  command time -f "CPU:%P  MEM(kbytes):%M" python3 DendropySingle.py gene_trees.1.$n.$k gene_trees.1.$n.$k "$m" &> gene_trees.1.$n.$k.ds.log
done
