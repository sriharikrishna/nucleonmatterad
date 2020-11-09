#make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 NUCMAT=1 CASE=pnm
cat pnm/av18uix_pnm.in >temp_$1_$2.dat
echo $1 >> temp_$1_$2_$3_$4.dat
echo $2 >> temp_$1_$2_$3_$4.dat
echo $3 >> temp_$1_$2_$3_$4.dat
echo $4 >> temp_$1_$2_$3_$4.dat
start=`date +%s`
time ./nmadv < temp_$1_$2_$3_$4.dat > out_tap_all_bfgs_pnm_$1_$2_$3_$4
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_bfgs_pnm_$1_$2_$3_$4
