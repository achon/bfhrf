n=48
for k in 1000 5000 10000 14446
do
  echo $k : $(date)
  command time -f "CPU:%P  MEM(kbytes):%M" python3 DendropySingle.py avian.48.$k avian.48.$k $n &> avian.48.$k.ds.log
done