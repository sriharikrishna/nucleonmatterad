#make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 CASE=snm CUSTOM_INPUTS=1
cat snm/av18uix.in >temp_$1_$2.dat
echo $1 $2 >> temp_$1_$2.dat 

start=`date +%s`
time ./nmadv < temp_$1_$2.dat > out_tap_all_dfo_snm_$1_$2
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_dfo_snm_$1_$2
