#make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 BFGS=1 CASE=snm
cat snm/av18uix.in circle.txt >temp_$1_$2.dat
echo $1 >> temp_$1_$2.dat
echo $2 >> temp_$1_$2.dat
echo $3 $4 >> temp_$1_$2.dat
echo $5 $6 >> temp_$1_$2.dat
start=`date +%s`
time ./nmadv < temp_$1_$2.dat > out_tap_all_bfgs_snm_$1_$2
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_bfgs_snm_$1_$2
