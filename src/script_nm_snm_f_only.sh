#mkdir -p snm/obj/; make clean; make -f MakefileTapf clean; make prep CASE=snm; make CASE=snm CUSTOM_INPUTS=1 NUCMAT=1
cat snm/av18uix.in >temp_$1_$2.dat
echo $1 $2 >> temp_$1_$2.dat 

sedify() {
  sed -i "s/rho=/rho=$1/g" $5
  sed -i "s/lc=/lc=$2/g" $5
  sed -i "s/ls=/ls=$3/g" $5
  sed -i "s/lt=/lt=$4/g" $5
}

sedify $3 $4 $5 $6  temp_$1_$2.dat

start=`date +%s`
time ./snm/snm.x < temp_$1_$2.dat > out_nm_snm_$1_$2
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_nm_snm_$1_$2
