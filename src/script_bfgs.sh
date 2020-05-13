make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=snm prep ; make -f MakefileTapf ALL=1 BFGS=1 CASE=snm
./nmadv < snm/av18uix.in > out_tap_all_bfgs
