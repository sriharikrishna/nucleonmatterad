
FILE=nmadv
if [ -f "$FILE" ]; then
    echo "$FILE exists. So not compiling"
else 
    echo "$FILE does not exist. Compiling"
    make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 NUCMAT=1 CASE=pnm VARS=9
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
time ./nmadv < temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$(printf %.3f $8)_$(printf %.3f $9)_${10}_${11}_${12}_${13}.dat > out_tap_all_nucmat_pnm_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$(printf %.3f $8)_$(printf %.3f $9)_${10}_${11}_${12}_${13}
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_nucmat_pnm_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$(printf %.3f $8)_$(printf %.3f $9)_${10}_${11}_${12}_${13}
