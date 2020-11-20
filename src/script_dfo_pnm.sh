#make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 CASE=pnm CUSTOM_INPUTS=1
cat pnm/av18uix_pnm.in >temp_$1_$2_$3_$4.dat
echo $1 $2 $3 $4 >> temp_$1_$2_$3_$4.dat 

start=`date +%s`
time ./nmadv < temp_$1_$2_$3_$4.dat > out_tap_all_dfo_pnm_$1_$2_$3_$4
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_dfo_pnm_$1_$2_$3_$4
