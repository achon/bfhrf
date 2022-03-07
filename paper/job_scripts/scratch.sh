# avian

n=48
cpu=2
for k in 1000 5000 10000 14446
do
  echo $k : $(date)
  command time -f "TIME:%E CPU:%P MEM(kb):%M" python3 DendropySingleMP.py avian.48.$k avian.48.$k $n $cpu
done