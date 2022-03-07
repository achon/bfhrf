#!/bin/bash

# path
base=/home/alvin/bfh/data

# commands
cd $base/avian
n=48
for k in 1000 5000 10000 14446
do
  echo avian $k : $(date)
  command time -f "TIME:%E CPU:%P MEM(kb):%M" ./hashrf avian.48.$k 0 -o avian.48.$k.hashrf
  echo avian $k : $(date)
done

cd $base/insect
for k in 1000 50000 100000 149278
do
  echo insect $k : $(date)
  command time -f "TIME:%E CPU:%P MEM(kb):%M" ./hashrf insect.144.$k 0 -o insect.144.$k.hashrf
  echo insect $k : $(date)
done

cd $base/k.99.100000
n=100
for k in 1000 25000 50000 75000 100000
do
  echo k $k : $(date)
  command time -f "TIME:%E CPU:%P MEM(kb):%M" ./hashrf gene_trees.1.99.$k 0 -o gene_trees.1.99.$k.hashrf
  echo k $k : $(date)
done

cd n $base/n.x.1000
k=1000
for n in 99 249 499 749 999
do
  ((m=n+1))
  echo "$m" : $(date)
  command time -f "TIME:%E CPU:%P MEM(kb):%M" ./hashrf gene_trees.1.$n.$k 0 -o gene_trees.1.$n.$k.hashrf
  echo k $k : $(date)
done
