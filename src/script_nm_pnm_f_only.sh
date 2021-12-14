mkdir -p pnm/obj/; make clean; make -f MakefileTapf clean; make prep CASE=pnm; make CASE=pnm CUSTOM_INPUTS=1 NUCMAT=1
cat pnm/av18uix_pnm.in >temp_$1_$2_$3_$4_$5_$6_$7_$8.dat
echo $1 $2 $3 $4 >> temp_$1_$2_$3_$4_$5_$6_$7_$8.dat 

sedify() {
  sed -i "s/ rho=/ rho=$1/g" $5
  sed -i "s/ lc=/ lc=$2/g" $5
  sed -i "s/ ls=/ ls=$3/g" $5
  sed -i "s/ lt=/ lt=$4/g" $5
}

sedify $5 $6 $7 $8 temp_$1_$2_$3_$4_$5_$6_$7_$8.dat

start=`date +%s`
time ./pnm/pnm.x < temp_$1_$2_$3_$4_$5_$6_$7_$8.dat > out_nm_pnm_$1_$2_$3_$4_$5_$6_$7_$8
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_nm_pnm_$1_$2_$3_$4_$5_$6_$7_$8
