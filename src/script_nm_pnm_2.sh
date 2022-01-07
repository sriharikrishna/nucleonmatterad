
FILE=nmadv
if [ -f "$FILE" ]; then
    echo "$FILE exists. So not compiling"
else 
    echo "$FILE does not exist. Compiling"
    make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 NUCMAT=1 CASE=pnm VARS=2
fi

cat pnm/search16.2.in >temp_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6.dat
echo $1 $2 $3 $4 >> temp_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6.dat

sedify() {
  sed -i "s/ rho=/ rho=$1/g" $5
  sed -i "s/ lc=/ lc=$2/g" $5
  sed -i "s/ ls=/ ls=$3/g" $5
  sed -i "s/ lt=/ lt=$4/g" $5
}

sedify $3 $4 $5 $6 temp_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6.dat

start=`date +%s`
time ./nmadv < temp_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6.dat > out_tap_all_nucmat_pnm_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_nucmat_pnm_$(printf %.3f $1)_$(printf %.3f $2)_$3_$4_$5_$6
