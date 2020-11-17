#mkdir -p pnm/obj/; make clean; make -f MakefileTapf clean; make prep CASE=pnm; make CASE=pnm CUSTOM_INPUTS=1
cat pnm/av18uix_pnm.in >temp.dat
echo $1 $2 $3 $4 >> temp.dat 

start=`date +%s`
time ./pnm/pnm.x < temp.dat > out_dfo_pnm_$1_$2_$3_4
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_dfo_pnm_$1_$2_$3_$4
