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
include Makefile-g95

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

OBJFIL= $(OBJ)/nmchain.o $(OBJ)/nmfts.o $(OBJ)/nmhnc.o $(OBJ)/nmmain.o \
   $(OBJ)/nmsub.o $(OBJ)/nmtbi.o

# parameter file names

PARAMS=nclude/params.$(CASE)

# subroutine object files (almost) independent of nucleus

OBJSUB= sub/dtime.o sub/headtime.o sub/linpack.o sub/minimi.o sub/numrec.o \
   sub/pot.o sub/eft_pot_r.o sub/dgamma.o

OBJ=$(CASE)/obj

all: $(CASE)/$(CASE).x

prep: 
	cp -p $(PARAMS) nclude/params.f

$(CASE)/$(CASE).x: $(CASE)/nm.o $(OBJFIL) $(OBJSUB)
	$(FC) $(OMP) $(CASE)/nm.o $(OBJFIL) $(OBJSUB) $(MLIB) -o $(CASE)/$(CASE).x

# case specific subroutines

$(CASE)/nm.o: nm.f $(PARAMS)
	$(FC) $(OMP) -c nm.f -o $(CASE)/nm.o

$(OBJ)/nmchain.o: nmchain.f $(PARAMS)
	$(FC) -c nmchain.f -o $(OBJ)/nmchain.o

$(OBJ)/nmfts.o: nmfts.f $(PARAMS)
	$(FC) -c nmfts.f -o $(OBJ)/nmfts.o

$(OBJ)/nmhnc.o: nmhnc.f $(PARAMS)
	$(FC) $(OMP) -c nmhnc.f -o $(OBJ)/nmhnc.o

$(OBJ)/nmmain.o: nmmain.f $(PARAMS)
	$(FC) -c nmmain.f -o $(OBJ)/nmmain.o

$(OBJ)/nmsub.o: nmsub.f $(PARAMS)
	$(FC) -c nmsub.f -o $(OBJ)/nmsub.o

$(OBJ)/nmtbi.o: nmtbi.f $(PARAMS)
	$(FC) -c nmtbi.f -o $(OBJ)/nmtbi.o

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
	rm -Rf sub/*.o pnm/*.o snm/*.o *.x *.o snm/*.x pnm/*.x *~ snm/*~ pnm/*~
