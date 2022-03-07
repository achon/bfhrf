n=144
for k in 1000 50000 100000 149278
do
  echo $k : $(date)
  command time -f "CPU:%P  MEM(kbytes):%M" python3 DendropySingle.py insect.144.$k insect.144.$k $n &> insect.144.$k.ds.log
done