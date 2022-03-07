# FROM ASTRAL FOR BHFS

# ml-genetrees.tre is the avian dataset.
# estimated_gene_trees.treesis the insect dataset.
# Renaming both for readability.
cp ml-genetrees.tre avian.48.14446
for k in 1000 5000 10000
do
  head -n $k avian.48.14446 > avian.48.$k
done

cp estimated_gene_trees.trees insect.144.149278
for k in 1000 50000 100000
do
  head -n $k insect.144.149278 > insect.144.$k
done


# Code to generate simphy data sets using ASTRAL parameters
# simphy has an extra species so n=10 is range 0:9 and input n=9
#k=1000, n=[100, 250, 500, 750, 1000]
rep=5
b=0.000001
t=2000000
genes=1000
for sp in 99 249 499 749 999
do
  ./simphy -rs $rep -rl U:$genes,$genes -rg 1 -st U:$t,$t -si U:1,1 -sl U:$sp,$sp -sb U:$b,$b -sp U:200000,200000 -hs LN:1.5,1 -hl LN:1.2,1 -hg LN:1.4,1 -su E:10000000 -so U:1,1 -v 1  -cs 293745 -o n.$sp.$genes
  rm n.$sp.$genes/*/[sl]_tree*.trees
  for r in `ls -d n.$sp.$genes/*`
  do
    sed -i "" -e "s/_0_0//g" $r/g_trees*.trees;
  done
  for r in `ls -d n.$sp.$genes/*`
  do
    cat $r/g_trees*.trees > $r/gene_trees.$sp.$genes;
    rm  $r/g_trees*.trees;
  done
  for i in 1 2 3 4 5
  do
    mv n.$sp.$genes/$i/gene_trees.$sp.$genes n.$sp.$genes/$i/gene_trees.$i.$sp.$genes
  done
done

#n=100, k=[1K, 25K, 50K, 75K, 100K]
#Generate genes=100000 then take subsets for smaller values of k
rep=5
b=0.000001
t=2000000
sp=19
genes=100000
./simphy -rs $rep -rl U:$genes,$genes -rg 1 -st U:$t,$t -si U:1,1 -sl U:$sp,$sp -sb U:$b,$b -sp U:200000,200000 -hs LN:1.5,1 -hl LN:1.2,1 -hg LN:1.4,1 -su E:10000000 -so U:1,1 -v 1  -cs 293745 -o k.$sp.$genes
rm k.$sp.$genes/*/[sl]_tree*.trees
for r in `ls -d k.$sp.$genes/*`
do
  echo $r/g_trees*.trees| xargs sed -i "" -e "s/_0_0//g";
done
for r in `ls -d k.$sp.$genes/*`
do
  echo $r
  echo $r/g_trees*.trees | xargs cat >> $r/gene_trees.$sp.$genes;
  echo $r/g_trees*.trees | xargs rm;
  for k in 1000 25000 50000 75000
  do
    head -n $k $r/gene_trees.$sp.$genes > $r/gene_trees.$sp.$k
    done
done
for i in 1 2 3 4 5
do
  for g in 1000 25000 50000 75000 100000
  do
    mv k.$sp.$genes/$i/gene_trees.$sp.$g k.$sp.$genes/$i/gene_trees.$i.$sp.$g
  done
done

