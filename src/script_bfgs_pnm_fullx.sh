#make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 BFGS=1 FULLX=1 CASE=pnm
cat pnm/av18uix_pnm.in circle.txt >temp_$1_$2.dat
echo $1 >> temp_$1_$2.dat
echo $2 >> temp_$1_$2.dat
echo $3 $4 >> temp_$1_$2.dat
echo $5 $6 >> temp_$1_$2.dat
echo $7 $8 >> temp_$1_$2.dat
echo $9 ${10} >> temp_$1_$2.dat
echo ${11} ${12} >> temp_$1_$2.dat
echo ${13} ${14} >> temp_$1_$2.dat
echo ${15} ${16} >> temp_$1_$2.dat
start=`date +%s`
time ./nmadv < temp_$1_$2.dat > out_tap_all_bfgs_pnm_$1_$2
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_bfgs_pnm_$1_$2
