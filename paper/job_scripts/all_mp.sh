#!/bin/bash

# path
base=/home/alvin/bfh/data

# commands
n=48
for cpu in 8 16
do
  cd $base/avian
  for k in 1000 5000 10000 14446
  do
    echo $k : $(date)
    command time -f "TIME:%E CPU:%P MEM(kb):%M" python3 DendropySingleMP.py avian.48.$k avian.48.$k $n $cpu &> avian.48.$k.$cpu.dsmp.log
    command time -f "TIME:%E CPU:%P MEM(kb):%M" python3 BipartitionFrequencyHash.py avian.48.$k avian.48.$k $n $cpu &> avian.48.$k.$cpu.bfh.log
  done

  cd $base/insect
  for k in 1000 50000 100000 149278
  do
    echo $k : $(date)
    command time -f "TIME:%E CPU:%P MEM(kb):%M" python3 DendropySingleMP.py insect.144.$k insect.144.$k $n $cpu &> insect.144.$k.$cpu.dsmp.log
    command time -f "TIME:%E CPU:%P MEM(kb):%M" python3 BipartitionFrequencyHash.py insect.144.$k insect.144.$k $n $cpu &> insect.144.$k.$cpu.bfh.log
  done

  cd $base/k.99.100000
  n=100
  for k in 1000 25000 50000 75000 100000
  do
    echo $k : $(date)
    command time -f "TIME:%E CPU:%P MEM(kb):%M" python3 DendropySingleMP.py gene_trees.1.99.$k gene_trees.1.99.$k $n $cpu &> gene_trees.1.99.$k.$cpu.dsmp.log
    command time -f "TIME:%E CPU:%P MEM(kb):%M" python3 BipartitionFrequencyHash.py gene_trees.1.99.$k gene_trees.1.99.$k $n $cpu &> gene_trees.1.99.$k.$cpu.bfh.log
  done

  cd $base/n.x.1000
  k=1000
  for n in 99 249 499 749 999
  do
    ((m=n+1))
    echo "$m" : $(date)
    command time -f "TIME:%E CPU:%P MEM(kb):%M" python3 DendropySingleMP.py gene_trees.1.$n.$k gene_trees.1.$n.$k "$m" $cpu &> gene_trees.1.$n.$k.$cpu.dsmp.log
    command time -f "TIME:%E CPU:%P MEM(kb):%M" python3 BipartitionFrequencyHash.py gene_trees.1.$n.$k gene_trees.1.$n.$k "$m" $cpu &> gene_trees.1.$n.$k.$cpu.bfh.log
  done
done