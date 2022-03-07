#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=128G
#SBATCH --time=2-0
#SBATCH --partition=bigram,biocrunch

# everything below this line is optional, but are nice to have quality of life things
#SBATCH --output=job.%J.out # tell it to store the output console text to a file called job.<assigned job number>.out
#SBATCH --error=job.%J.err # tell it to store the error messages from the program (if it doesn't write them to normal console output) to a file called job.<assigned job muber>.err
#SBATCH --job-name="run_bfh_n_48" # a nice readable name to give your job so you know what it is when you see it in the queue, instead of just numbers

# modules

# path

cd /work/LAS/xqhuang-lab/achon/bfh
cd data/n.x.1000

# commands
source /work/LAS/xqhuang-lab/achon/p3-venv/bin/activate #enter venv

k=1000
cpu=48
for n in 99 249 499 749 999
do
  ((m=n+1))
  echo "$m" : $(date)
  command time -f "CPU:%P  MEM(kbytes):%M" python3 BipartitionFrequencyHash.py gene_trees.1.$n.$k gene_trees.1.$n.$k "$m" $cpu &> gene_trees.1.$n.$k.bfh.$cpu.log
done

