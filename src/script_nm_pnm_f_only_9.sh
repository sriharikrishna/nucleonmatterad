FILE=./pnm/pnm.x
if [ -f "$FILE" ]; then
    echo "$FILE exists. So not compiling"
else 
    echo "$FILE does not exist. Compiling"
    mkdir -p pnm/obj/; make clean; make -f MakefileTapf clean; make prep CASE=pnm; make CASE=pnm CUSTOM_INPUTS=1 
fi
cat pnm/search16.9.in >temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$(printf %.3f $8)_$(printf %.3f $9)_${10}_${11}_${12}_${13}.dat
echo $1 $2 $3 $4 $5 $6 $7 $8 $9 >> temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$(printf %.3f $8)_$(printf %.3f $9)_${10}_${11}_${12}_${13}.dat 

sedify() {
  sed -i "s/ rho=/ rho=$1/g" $5
  sed -i "s/ lc=/ lc=$2/g" $5
  sed -i "s/ ls=/ ls=$3/g" $5
  sed -i "s/ lt=/ lt=$4/g" $5
}

sedify ${10} ${11} ${12} ${13} temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$(printf %.3f $8)_$(printf %.3f $9)_${10}_${11}_${12}_${13}.dat

start=`date +%s`
time ./pnm/pnm.x < temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$(printf %.3f $8)_$(printf %.3f $9)_${10}_${11}_${12}_${13}.dat > out_nm_pnm_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$(printf %.3f $8)_$(printf %.3f $9)_${10}_${11}_${12}_${13}
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_nm_pnm_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$(printf %.3f $8)_$(printf %.3f $9)_${10}_${11}_${12}_${13}
