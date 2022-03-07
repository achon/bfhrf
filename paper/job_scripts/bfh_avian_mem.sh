#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=32G
#SBATCH --time=2-0
#SBATCH --partition=biocrunch

# everything below this line is optional, but are nice to have quality of life things
#SBATCH --output=job.%J.out # tell it to store the output console text to a file called job.<assigned job number>.out
#SBATCH --error=job.%J.err # tell it to store the error messages from the program (if it doesn't write them to normal console output) to a file called job.<assigned job muber>.err
#SBATCH --job-name="m_avian" # a nice readable name to give your job so you know what it is when you see it in the queue, instead of just numbers

# modules

# path

cd /work/LAS/xqhuang-lab/achon/bfh
cd data/avian

# commands
source /work/LAS/xqhuang-lab/achon/p3-venv/bin/activate #enter venv
n=48
for cpu in 24 48
do
  for k in 1000 5000 10000 14446
  do
    echo $k : $(date)
    command time -f "CPU:%P  MEM(kbytes):%M" python3 BipartitionFrequencyHash.py avian.48.$k avian.48.$k $n $cpu &> avian.48.$k.bfh.$cpu.mem.log
  done
done