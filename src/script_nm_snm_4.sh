
FILE=nmadv
if [ -f "$FILE" ]; then
    echo "$FILE exists. So not compiling"
else 
    echo "$FILE does not exist. Compiling"
    make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 NUCMAT=1 CASE=snm VARS=4
fi

cat snm/search.in20 >temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$5_$6_$7_$8.dat
echo $1 $2 $3 $4 >> temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$5_$6_$7_$8.dat

sedify() {
  sed -i "s/ rho=/ rho=$1/g" $5
  sed -i "s/ lc=/ lc=$2/g" $5
  sed -i "s/ ls=/ ls=$3/g" $5
  sed -i "s/ ll=/ ll=$4/g" $5
}

sedify $5 $6 $7 $8 temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$5_$6_$7_$8.dat

start=`date +%s`
time ./nmadv < temp_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$5_$6_$7_$8.dat > out_tap_all_nucmat_snm_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$5_$6_$7_$8
end=`date +%s`
runtime=$((end-start))

echo runtime $runtime >> out_tap_all_nucmat_snm_$(printf %.3f $1)_$(printf %.3f $2)_$(printf %.3f $3)_$(printf %.3f $4)_$5_$6_$7_$8
