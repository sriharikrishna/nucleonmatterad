make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf DIFF=-DDO_DOR CASE=pnm
./nmad < pnm/av18uix_pnm.in > out_tap_pnm_1_bls0_dor
make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf DIFF=-DDO_BST CASE=pnm
./nmad < pnm/av18uix_pnm.in > out_tap_pnm_1_bls0_bst
make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf DIFF=-DDO_BTN CASE=pnm
./nmad < pnm/av18uix_pnm.in > out_tap_pnm_1_bls0_btn
make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf DIFF=-DDO_BLS CASE=pnm
./nmad < pnm/av18uix_pnm.in > out_tap_pnm_1_bls0_bls
make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf DIFF=-DDO_AST CASE=pnm
./nmad < pnm/av18uix_pnm.in > out_tap_pnm_1_bls0_ast
make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf DIFF=-DDO_ATN CASE=pnm
./nmad < pnm/av18uix_pnm.in > out_tap_pnm_1_bls0_atn
make -f MakefileTapf clean; make -f Makefile clean; make -f Makefile CASE=pnm prep ; make -f MakefileTapf DIFF=-DDO_ALS CASE=pnm
./nmad < pnm/av18uix_pnm.in > out_tap_pnm_1_bls0_als

