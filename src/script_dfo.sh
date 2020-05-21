
#make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 CASE=snm numCore_dv.pre.f
#exit
make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 CASE=snm
cat snm/av18uix.in circle.txt >temp.dat
echo $1 >> temp.dat 
echo $2 >> temp.dat 
./nmadv < temp.dat > out_tap_all_dfo
