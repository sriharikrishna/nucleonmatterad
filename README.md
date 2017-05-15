# nucleonmatterad
Algorithmic differentiation of either symmetric nuclean matter (SNM) or putre nucleon matter (PNM). 

file		subroutines	calls           notes
                (entries)
----		-----------	-----           -----
nm.f		nm		header          driving code
		nucmat		timer           computes energy for one set of
		(nminit)	minimi          parameters and optionally 
		(nmfin)		nmmain          calls minimization package
				nmout           to search for best values

nmspot.f        nmspot          header          alternate driving code 
                                timer		for single-particle potential
                                nmmain

nm.in                                           sample input deck

nm.out                                          sample output deck

nmmain.f:	nmmain		nmfts           main energy computation
				nmhot           calls correlation generator
				nmsps           and fhnc/soc chain solver
                                nmchain         sums 2-body and separable
				nmhnc           energies
				nmtbi
				nmpion

nmchain.f       nmchain		-----           sums chain energies

nmfts.f		nmfts		setpot          correlation generator
				nmhot
				pot (4x)

nmhnc.f		nmhnc		locate          solves fhnc/soc equations
                                polint

nmtbi.f		nmtbi		-----           Vijk, U, UF integrations

nmsub.f		nmhot		-----           finite temperature setup
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
