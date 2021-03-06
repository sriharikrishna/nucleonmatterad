# Makefile for nuclear/neutron matter
#
# Invoke this file with the commands:
# for symmetric nuclear matter (snm)
# make (-n) prep CASE=snm : sets up nclude files
# make (-n) CASE=snm      : compiles for snm.x
# or for pure neutron matter (pnm)
# make (-n) prep CASE=pnm : sets up nclude files
# make (-n) CASE=pnm      : compiles for pnm.x
#

#include Makefile-theory
#include Makefile-g95
include Makefile-gfortran

#
# may invoke OPENMP=on
#

ifeq ($(OPENMP), on)
   OMP=$(OMPFLAGS)
else
   OMP=
endif

# required object files
# you must define CASE either as environment variables or on
# the command line

OBJ=$(CASE)/obj

OBJFIL= $(OBJ)/nmchain.o $(OBJ)/nmfts.o $(OBJ)/nmhnc.o $(OBJ)/nmmain.o \
   $(OBJ)/nmsub.o $(OBJ)/nmtbi.o $(OBJ)/params.o $(OBJ)/nmsubmod.o $(OBJ)/nmvar.o $(OBJ)/nmvarcopypassive.o 

# parameter file names

PARAMS=nclude/params.$(CASE)

# subroutine object files (almost) independent of nucleus

OBJSUB= sub/dtime.o sub/headtime.o sub/linpack.o sub/minimi.o sub/numrec.o \
   sub/pot.o sub/eft_pot_r.o sub/dgamma.o 


all: $(CASE)/$(CASE).x

prep: 
	cp -p $(PARAMS) nclude/params.f

$(OBJ)/params.o: nclude/params.f
	$(FC) $(OMP) -c nclude/params.f -o $(OBJ)/params.o

$(CASE)/$(CASE).x: $(CASE)/nm.o $(OBJFIL) $(OBJSUB)
	$(FC) $(OMP) $(CASE)/nm.o $(OBJFIL) $(OBJSUB) $(MLIB) -o $(CASE)/$(CASE).x

# case specific subroutines

$(CASE)/nm.o: nm.f $(OBJ)/params.o $(OBJ)/nmvar.o $(OBJ)/nmvarcopypassive.o
	$(FC) $(OMP) -c nm.f -o $(CASE)/nm.o

$(OBJ)/nmchain.o: nmchain.f $(OBJ)/params.o $(OBJ)/nmsubmod.o
	$(FC) -c nmchain.f -o $(OBJ)/nmchain.o

$(OBJ)/nmfts.o: nmfts.f $(OBJ)/params.o
	$(FC) -c nmfts.f -o $(OBJ)/nmfts.o

$(OBJ)/nmhnc.o: nmhnc.f $(OBJ)/params.o
	$(FC) $(OMP) -c nmhnc.f -o $(OBJ)/nmhnc.o

$(OBJ)/nmmain.o: nmmain.f $(OBJ)/params.o $(OBJ)/nmhnc.o $(OBJ)/nmtbi.o $(OBJ)/nmsubmod.o
	$(FC) -c nmmain.f -o $(OBJ)/nmmain.o

$(OBJ)/nmsub.o: nmsub.f $(OBJ)/params.o
	$(FC) -c nmsub.f -o $(OBJ)/nmsub.o

$(OBJ)/nmsubmod.o: nmsubmod.f $(OBJ)/params.o
	$(FC) -c nmsubmod.f -o $(OBJ)/nmsubmod.o

$(OBJ)/nmtbi.o: nmtbi.f $(OBJ)/params.o
	$(FC) -c nmtbi.f -o $(OBJ)/nmtbi.o

$(OBJ)/nmvar.o: nmvar.f $(OBJ)/params.o
	$(FC) -c nmvar.f -o $(OBJ)/nmvar.o

$(OBJ)/nmvarcopypassive.o: nmvar.f $(OBJ)/params.o $(OBJ)/nmvar.o
	$(FC) -c nmvarcopypassive.f -o $(OBJ)/nmvarcopypassive.o

# common subroutines

sub/dtime.o: sub/dtime.c
	gcc -c sub/dtime.c -o sub/dtime.o

sub/headtime.o: sub/headtime.f
	$(FC) -c sub/headtime.f -o sub/headtime.o

sub/linpack.o: sub/linpack.f
	$(FC) -c sub/linpack.f -o sub/linpack.o

sub/minimi.o: sub/minimi.f
	$(FC) -c sub/minimi.f -o sub/minimi.o

sub/numrec.o: sub/numrec.f
	$(FC) -c sub/numrec.f -o sub/numrec.o

sub/pot.o: sub/pot.f
	$(FC) -c sub/pot.f -o sub/pot.o

sub/eft_pot_r.o: sub/eft_pot_r.f sub/params_new_pot.f
	$(FC) -c sub/eft_pot_r.f -o sub/eft_pot_r.o

sub/dgamma.o: sub/dgamma.f
	$(FC) -c sub/dgamma.f -o sub/dgamma.o

clean:
	rm -Rf $(OBJFIL) $(OBJSUB) sub/*.o pnm/*.o snm/*.o *.x *.o snm/*.x pnm/*.x *~ snm/*~ pnm/*~ *.mod snm/obj/*.mod snm/obj/*.o pnm/obj/*.mod pnm/obj/*.o
