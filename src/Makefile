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
   $(OBJ)/nmsub.o $(OBJ)/nmtbi.o $(OBJ)/isnan_d.o 

# parameter file names

PARAMS=nclude/params.$(CASE)

# subroutine object files (almost) independent of nucleus

OBJSUB= sub/dtime.o sub/headtime.o sub/linpack.o sub/minimi.o sub/numrec.pre.o \
   sub/pot.o sub/eft_pot_r.pre.o sub/dgamma.pre.o 


all: $(CASE)/$(CASE).x

prep: 
	cp -p $(PARAMS) nclude/params.f


$(CASE)/$(CASE).x: $(CASE)/nm.pre.o $(OBJFIL) $(OBJSUB)
	$(FC) $(OMP) /pwd/$(CASE)/nm.pre.o $(addprefix /pwd/,$(OBJFIL)) $(addprefix /pwd/,$(OBJSUB)) $(MLIB) -o /pwd/$(CASE)/$(CASE).x

$(CASE)/nm.pre.o: nm.pre.f 
	$(FC)  -c /pwd/nm.pre.f -o /pwd/$(CASE)/nm.pre.o

# case specific subroutines

$(OBJ)/nmchain.o: nmchain.f 
	$(FC) -c /pwd/nmchain.f -o /pwd/$(OBJ)/nmchain.o

$(OBJ)/nmfts.o: nmfts.f 
	$(FC) -c /pwd/nmfts.f -o /pwd/$(OBJ)/nmfts.o

$(OBJ)/nmhnc.o: nmhnc.f 
	$(FC) $(OMP) -c /pwd/nmhnc.f -o /pwd/$(OBJ)/nmhnc.o

$(OBJ)/nmmain.o: nmmain.f $(OBJ)/nmhnc.o $(OBJ)/nmtbi.o 
	$(FC) -c /pwd/nmmain.f -o /pwd/$(OBJ)/nmmain.o

$(OBJ)/nmsub.o: nmsub.f 
	$(FC) -c /pwd/nmsub.f -o /pwd/$(OBJ)/nmsub.o

$(OBJ)/nmtbi.o: nmtbi.f 
	$(FC) -c /pwd/nmtbi.f -o /pwd/$(OBJ)/nmtbi.o

$(OBJ)/nmvar.o: nmvar.f 
	$(FC) -c /pwd/nmvar.f -o /pwd/$(OBJ)/nmvar.o

$(OBJ)/isnan_d.o: isnan_d.f 
	$(FC) -c /pwd/isnan_d.f -o /pwd/$(OBJ)/isnan_d.o

# common subroutines

sub/dtime.o: sub/dtime.c
	echo "hello"
	echo docker run -v `pwd`:/pwd wsmoses/flang clang -flto -fuse-ld=lld-7 -c /pwd/sub/dtime.c -o /pwd/sub/dtime.o
	docker run -v `pwd`:/pwd wsmoses/flang clang -flto -fuse-ld=lld-7 -c /pwd/sub/dtime.c -o /pwd/sub/dtime.o

sub/headtime.o: sub/headtime.f
	$(FC) -c /pwd/sub/headtime.f -o /pwd/sub/headtime.o

sub/linpack.o: sub/linpack.f
	$(FC) -c /pwd/sub/linpack.f -o /pwd/sub/linpack.o

sub/minimi.o: sub/minimi.f
	$(FC) -c /pwd/sub/minimi.f -o /pwd/sub/minimi.o

sub/numrec.o: sub/numrec.f
	$(FC) -c /pwd/sub/numrec.f -o /pwd/sub/numrec.o

sub/pot.o: sub/pot.f
	$(FC) -c /pwd/sub/pot.f -o /pwd/sub/pot.o

sub/eft_pot_r.o: sub/eft_pot_r.f /pwd/sub/params_new_pot.f
	$(FC) -c /pwd/sub/eft_pot_r.f -o /pwd/sub/eft_pot_r.o

CPPFLAGS= -DALLOW_DATA 
ifdef CUSTOM_INPUTS
CPPFLAGS+=-DCUSTOM_INPUTS
endif
ifdef NUCMAT
CPPFLAGS+=-DONLY_NUCMAT
endif
ifeq ($(CASE), snm)
  CPPFLAGS+=-DCASE_SNM
else
  CPPFLAGS+=-DCASE_PNM
endif
CPPCMD = cat $< |  cpp $(CPPFLAGS) -I../ -I./ -traditional-cpp -P 

%.pre.f: %.f
	$(CPPCMD) > $@


sub/dgamma.o: sub/dgamma.f
	$(FC) -c /pwd/sub/dgamma.f -o /pwd/sub/dgamma.o

clean:
	rm -Rf $(OBJFIL) $(OBJSUB) sub/*.o pnm/*.o snm/*.o *.x *.o snm/*.x pnm/*.x *~ snm/*~ pnm/*~ *.mod snm/obj/*.mod snm/obj/*.o pnm/obj/*.mod pnm/obj/*.o
