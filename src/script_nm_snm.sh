mkdir -p snm/obj/; make clean; make -f MakefileTapf clean; make prep CASE=snm; make CASE=snm CUSTOM_INPUTS=1 NUCMAT=1
cat snm/av18uix.in >temp.dat
echo $1 $2 >> temp.dat 

start=`date +%s`
time ./snm/snm.x < temp.dat > out_nm_snm_$1_$2
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_nm_snm_$1_$2
