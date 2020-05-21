
#make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 CASE=snm numCore_dv.pre.f
#exit
make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf ALL=1 CASE=pnm
./nmadv < pnm/av18uix_pnm.in > out_tap_all_dfo_pnm
