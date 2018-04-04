# nucleonmatterad
Algorithmic differentiation of either symmetric nuclean matter (SNM) or pure nucleon matter (PNM). 

Description of the files in the original code:

file		subroutines	calls           notes
                (entries)
----		-----------	-----           -----
src/nm.f	nm		header          driving code
		nucmat		timer           computes energy for one set of
		(nminit)	minimi          parameters and optionally 
		(nmfin)		nmmain          calls minimization package
				nmout           to search for best values

src/nmspot.f    nmspot          header          alternate driving code 
                                timer		for single-particle potential
                                nmmain

nm.in                                           sample input deck

nm.out                                          sample output deck

src/nmmain.f:	nmmain		nmfts           main energy computation
				nmhot           calls correlation generator
				nmsps           and fhnc/soc chain solver
                                nmchain         sums 2-body and separable
				nmhnc           energies
				nmtbi
				nmpion

src/nmchain.f	nmchain		-----           sums chain energies

src/nmfts.f	nmfts		setpot          correlation generator
				nmhot
				pot (4x)

src/nmhnc.f	nmhnc		locate          solves fhnc/soc equations
                                polint

src/nmtbi.f	nmtbi		-----           Vijk, U, UF integrations

src/nmsub.f	nmhot		-----           finite temperature setup
		nmsps                           single-particle potential setup
		nmpion                          excess pion calculation
		nmout                           output routine
		ac                              operator matrix functions
		(acex)
		(acl2)
		(acl2ex)
		(al2)

pot.f		setpot		-----           potential initialization
		(pot)                           potential subroutine

headtime.f      header                          output header  | machine- 
		timer                           timing routine | dependent

minimi.f	minimi		dsvdc           minimization package for
		monitr		dgeco           finding optimal parameters
		functn		dgesl
		smplex                          simplex routine
		versrt
		centrd
		reflct
		expand
		cntrct
		movevt
		checkt
		newset
		quadr                           quadratic routine
		indexx

numrec.f        locate          -----           numerical recipes
                polint

linpack.f	dgeco		                standard linpack routines
		dgedi                           needed by minimization
		dgefa                           package (use single-precision
		dsvdc                           variant on a cray)
		dasum
		daxpy
		ddot
		dnrm2
		drot
		drotg
		dscal
		dswap
		idamax



Obtaining the Nucleonmatter code
--------------------------------
The code can be cloned from
git clone https://github.com/sriharikrishna/nucleonmatterad.git

or checked out of
svn co https://repocafe.cels.anl.gov/repos/nucleonmatter

Installing OpenAD
-----------------
Comprehensive instructions on obtianing and building
OpenAD are available at
http://www.mcs.anl.gov/OpenAD/


Creating the differentiated code
-----------------
1. Every session must start with executing
source setenv.sh
OR 
source setenv.csh 
within OpenAD source directory. 

2. Change directory to src

3. Based on whether one wants to run pnm or snm versions and
   forward mode (tangent linear) or reverse mode (adjoint) AD 
   one can employ the following commands

   a.  Forward mode AD code for snm
   make -f MakefileOpenADF clean ; make prep CASE=snm ; make -f MakefileOpenADF

   b.  Reverse mode AD code for snm
   make -f MakefileOpenAD clean ; make prep CASE=snm ; make -f MakefileOpenAD

   c.  Forward mode AD code for pnm
   make -f MakefileOpenADF clean ; make prep CASE=pnm ; make -f MakefileOpenADF

   d.  Reverse mode AD code for pnm
   make -f MakefileOpenADclean ; make prep CASE=pnm ; make -f MakefileOpenAD

4. Running the differentiated code
   -----------------
   Based on the differenitated code, one should either
   nmad a < inputfile > outputfile
   or
   nmad a < inputfile > outputfile

Running Finite Differences
----------------- 
1. Change directory to src
2. Run either
   make clean ; make prep CASE=snm ; make CASE=snm
   or
   make clean ; make prep CASE=pnm ; make CASE=pnm
3. Then run
   snm/snm.x f < inputfile > outputfile
   or
   pnm/pnm.x f < inputfile > outputfile
