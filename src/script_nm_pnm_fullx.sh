make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 FULLX=1 NUCMAT=1 CASE=pnm
cat pnm/av18uix_pnm.in >temp_$1_$2_$3_$4_$5.dat
echo $1 $2 $3 $4 $5 >> temp_$1_$2_$3_$4_$5.dat
start=`date +%s`
time ./nmadv < temp_$1_$2_$3_$4_$5.dat > out_tap_all_nucmat_pnm_$1_$2_$3_$4_$5
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_nucmat_pnm_$1_$2_$3_$4_$5
