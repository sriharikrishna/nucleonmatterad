#mkdir -p pnm/obj/; make clean; make -f MakefileTapf clean; make prep CASE=pnm; make CASE=pnm CUSTOM_INPUTS=1 NUCMAT=1
cat pnm/av18uix_pnm.in >temp_$1_$2_$3_$4.dat
echo $1 $2 $3 $4 >> temp_$1_$2_$3_$4.dat 

start=`date +%s`
time ./pnm/pnm.x < temp_$1_$2_$3_$4.dat > out_nm_pnm_$1_$2_$3_$4
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_nm_pnm_$1_$2_$3_$4
