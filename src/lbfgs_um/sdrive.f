C
C     ***********************
C     SIMPLE DRIVER FOR LBFGS
C     ***********************
C
C     Example of driver for LBFGS routine, using a
C     simple test problem. The solution point is at 
C     X=(1,...,1) and the optimal function value of 0.
C
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
#ifndef DO_FULLX
      SUBROUTINE SDRIVE(NDIM,MSAVE,X,funk)
#else
      SUBROUTINE SDRIVE(NDIM,NCASE,MSAVE,X,funk)
#endif
      INTEGER NDIM,MSAVE,NCASE
      parameter (nlog=0,nin=5,nout=6)
C      PARAMETER(NDIM=2000,MSAVE=7,NWORK=NDIM*(2*MSAVE +1)+2*MSAVE)
      DOUBLE PRECISION X(NDIM),G(NDIM),DIAG(NDIM)
      DOUBLE PRECISION W((NDIM)*(2*MSAVE +1)+2*MSAVE)
      DOUBLE PRECISION F,EPS,XTOL,GTOL,T1,T2,STPMIN,STPMAX
      INTEGER IPRINT(2),IFLAG,ICALL,N,M,MP,LP,J
      LOGICAL DIAGCO
      DOUBLE PRECISION FLOCALD
      logical lprt
C
C     The driver for LBFGS must always declare LB2 as EXTERNAL
C
      EXTERNAL LB2
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
C
      N=NDIM
      M=MSAVE
      IPRINT(1)= 1
      IPRINT(2)= 3
C
C     We do not wish to provide the diagonal matrices Hk0, and 
C     therefore set DIAGCO to FALSE.
C
      DIAGCO= .FALSE.
c      EPS= 1.0D-5
c      EPS= 1.0D-4
c      EPS= 1.0D-2
      EPS= 0.020
      XTOL= 1.0D-16
      ICALL=0
      IFLAG=0
C
      
 20   CONTINUE
#ifndef DO_FULLX
      do i=1,ndim
        write(nlog,*) "x before flocald%d", G(i)
        write(nout,*) "x before flocald%d", G(i)
      end do
      G(1:ndim) =0.0
      call funk(x,NDIM,F,G,ICALL+1)
#else
      call funk(x,NCASE,F,G,ICALL+1)
#endif
      write(nlog,*) "NDIM", NDIM, "IFLAG", IFLAG
      write(nout,*) "NDIM", NDIM, "IFLAG", IFLAG
      write(nlog,*) "sdrive F", F
      write(nout,*) "sdrive F", F
      do i=1,ndim
        write(nlog,*) "sdrive G", G(i)
        write(nout,*) "sdrive G", G(i)
      end do
      if (isnan(F)) then
        stop "Stop NaN"
      end if
      CALL LBFGS(NDIM,M,X,F,G,DIAGCO,DIAG,IPRINT,
     &            EPS,XTOL,W,IFLAG)

      write(nlog,*) "IFLAG", IFLAG
      write(nout,*) "IFLAG", IFLAG
      IF(IFLAG.LE.0) GO TO 50
      ICALL=ICALL + 1
C     We allow at most 2000 evaluations of F and G
      IF(ICALL.GT.2000) GO TO 50
      write(nlog,9997) ICALL,F,(x(i),i=1,ndim)
      write(nout,9997) ICALL,F,(x(i),i=1,ndim)
 9997 format (' call # ',i3,' gives',f30.17,' at',(10f30.17))
      
      GO TO 20
  50  CONTINUE
      write(nlog,8128) icall,F,(x(i),i=1,ndim)
      write(nout,8128) icall,F,(x(i),i=1,ndim)
 8128 format ('after ',i3,' function calls final value is',
     & f12.5,' at',(10f30.17) )
      END SUBROUTINE
C
C     ** LAST LINE OF SIMPLE DRIVER (SDRIVE) **

