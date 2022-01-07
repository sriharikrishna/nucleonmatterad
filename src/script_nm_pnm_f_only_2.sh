FILE=./pnm/pnm.x
if [ -f "$FILE" ]; then
    echo "$FILE exists. So not compiling"
else 
    echo "$FILE does not exist. Compiling"
    mkdir -p pnm/obj/; make clean; make -f MakefileTapf clean; make prep CASE=pnm; make CASE=pnm CUSTOM_INPUTS=1 
fi
cat pnm/search16.2.in >temp_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6.dat
echo $1 $2 >> temp_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6.dat 

sedify() {
  sed -i "s/ rho=/ rho=$1/g" $5
  sed -i "s/ lc=/ lc=$2/g" $5
  sed -i "s/ ls=/ ls=$3/g" $5
  sed -i "s/ lt=/ lt=$4/g" $5
}

sedify $3 $4 $5 $6 temp_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6.dat

start=`date +%s`
time ./pnm/pnm.x < temp_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6.dat > out_nm_pnm_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_nm_pnm_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6
