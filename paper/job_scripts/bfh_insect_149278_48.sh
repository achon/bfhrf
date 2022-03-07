#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=128G
#SBATCH --time=3-0
#SBATCH --partition=bigram,biocrunch

# everything below this line is optional, but are nice to have quality of life things
#SBATCH --output=job.%J.out # tell it to store the output console text to a file called job.<assigned job number>.out
#SBATCH --error=job.%J.err # tell it to store the error messages from the program (if it doesn't write them to normal console output) to a file called job.<assigned job muber>.err
#SBATCH --job-name="run_bfh_insect" # a nice readable name to give your job so you know what it is when you see it in the queue, instead of just numbers

# modules

# path

cd /work/LAS/xqhuang-lab/achon/bfh
cd data/insect

# commands
source /work/LAS/xqhuang-lab/achon/p3-venv/bin/activate #enter venv

n=144
cpu=48
k=149278
echo $k : $(date)
command time -f "CPU:%P  MEM(kbytes):%M" python3 BipartitionFrequencyHash.py insect.144.$k insect.144.$k $n $cpu &> insect.144.$k.bfh.$cpu.log
