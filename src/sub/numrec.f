c *id* locate **********************************************************
c given an array xx of length n, and given a value of x, reutrns a value
c of j such that x is between xx(j) and xx(j+1).  xx must be monotonic,
c either increasing or decreasing.  j=0 or j=n is returned to indicate
c that x is out of range
c ----------------------------------------------------------------------
      subroutine locate(xx,n,x,j)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension xx(n)
      jl=0
      ju=n+1
   10 if (ju-jl.gt.1) then
        jm=(ju+jl)/2
        if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
          jl=jm
        else
          ju=jm
        end if
      go to 10
      end if
      j=jl
      return
      end
c *id* polint **********************************************************
c given arrays xa and ya, each of length n, and given a value x, this
c routine returns a value y, and an error estimate dy.  if p(x) is the
c polynomial of degree n-1 such that p(xa(i)) = ya(i), i=1,...,n, then
c the returned value y = p(x).
c ----------------------------------------------------------------------
      subroutine polint(xa,ya,n,x,y,dy)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (nmax=10)
      dimension xa(n),ya(n),c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        end if
        c(i)=ya(i)
        d(i)=ya(i)
   11 continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if (den.eq.0.) pause
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
   12   continue
        if (2*ns.lt.n-m) then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        end if
        y=y+dy
   13 continue
      return
      end
