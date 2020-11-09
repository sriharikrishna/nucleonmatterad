make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 FULLX=1 NUCMAT=1 CASE=snm
cat snm/av18uix.in >temp_$1_$2_$3_$4_$5_$6_$7.dat
echo $1 $2 $3 $4 $5 $6 $7 >> temp_$1_$2_$3_$4_$5_$6_$7.dat
start=`date +%s`
time ./nmadv < temp_$1_$2_$3_$4_$5_$6_$7.dat > out_tap_all_nucmat_snm_$1_$2_$3_$4_$5_$6_$7
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_nucmat_snm_$1_$2_$3_$4_$5_$6_$7
