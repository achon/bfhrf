#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=16G
#SBATCH --time=0-6
#SBATCH --partition=bigram

# everything below this line is optional, but are nice to have quality of life things
#SBATCH --output=job.%J.out # tell it to store the output console text to a file called job.<assigned job number>.out
#SBATCH --error=job.%J.err # tell it to store the error messages from the program (if it doesn't write them to normal console output) to a file called job.<assigned job muber>.err
#SBATCH --job-name="run_bfh_insect" # a nice readable name to give your job so you know what it is when you see it in the queue, instead of just numbers

# modules

# path

cd /work/LAS/xqhuang-lab/achon/bfh
cd data/k.99.100000

# commands
source /work/LAS/xqhuang-lab/achon/p3-venv/bin/activate #enter venv

n=100
cpu=48
k=1000
echo $k : $(date)
command time -f "CPU:%P  MEM(kbytes):%M" python3 BipartitionFrequencyHash.py gene_trees.1.99.$k gene_trees.1.99.$k $n $cpu &> gene_trees.1.99.$k.bfh.$cpu.log

