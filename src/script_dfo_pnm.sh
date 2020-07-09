make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 CASE=pnm
cat pnm/av18uix_pnm.in circle.txt >temp.dat
echo $1 >> temp.dat 
echo $2 >> temp.dat 

start=`date +%s`
time ./nmadv < temp.dat > out_tap_all_dfo_pnm_$1_$2
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_dfo_pnm_$1_$2
