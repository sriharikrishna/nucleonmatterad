#make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 CASE=snm
cat snm/av18uix.in circle.txt >temp.dat
echo $1 >> temp.dat 
echo $2 >> temp.dat 

start=`date +%s`
time ./nmadv < temp.dat > out_tap_all_dfo_snm_$1_$2
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_dfo_snm_$1_$2
