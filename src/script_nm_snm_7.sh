
FILE=nmadv
if [ -f "$FILE" ]; then
    echo "$FILE exists. So not compiling"
else 
    echo "$FILE does not exist. Compiling"
    make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 NUCMAT=1 CASE=snm VARS=7
fi

cat snm/av6p.seven08.in  >temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$8_$9_${10}_${11}.dat
echo $1 $2 $3 $4 $5 $6 $7 $8 $9 >> temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$8_$9_${10}_${11}.dat 

sedify() {
  sed -i "s/ rho=/ rho=$1/g" $5
  sed -i "s/ lc=/ lc=$2/g" $5
  sed -i "s/ ls=/ ls=$3/g" $5
  sed -i "s/ lt=/ lt=$4/g" $5
}

sedify ${8} ${9} ${10} ${11} temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$8_$9_${10}_${11}.dat

start=`date +%s`
time ./nmadv < temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$8_$9_${10}_${11}.dat > out_tap_all_nucmat_snm_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$8_$9_${10}_${11}
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_nucmat_snm_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$(printf %.3f $5)_$(printf %.3f $6)_$(printf %.3f $7)_$8_$9_${10}_${11}
