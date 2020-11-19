c *id* minimi **********************************************************
c minimization package
c
c     find the unconstrained minimum of a function using the simplex
c     method followed by n-dimension quadratic fits
c     see J. A. Nelder & R. Mead, Computer Journal 7 (1965) 308.
c
c     maxcl1 - max number of simplex function values to use
c     tol1 - tolerance to obtain by simplex method.  the rms
c            deviation of the values in the current simplex must
c            be less than tol1 to stop.  note tol1 is absolute.
c     maxcl2 - max number of quadratic function values to use
c     tol2 - tolerance to obtain by quadratic method.  the absolute
c            value of the difference of the current quadratic
c            function value and the best value must be less than
c            tol2 to stop.  tol2 is absolute.
c     a - coefficent for reflection - suggested value is 1.
c     b - coefficent for contraction - suggested value is .5
c     c - coefficent for expansion -  suggested value is 2.
c     scale - an array of dimension n containing scale factors
c             such that  x(i)*scale(i) ^ 10.
c     x - an array of dimension n containing:
c         on input the initial values of the search parameters
c         on output the final values of the parameters.
c     fbest - will be set to the function value at x.
c     funk - the name of a subroutine that will evaluate the function
c            to be minimized as follows:
c                 subroutine funk ( x, n, f )
c            in which
c                 x is the position
c                 n is the number of parameters in x
c                 f must be set to    f(x)
c            funk must be declaired external in the calling program
c     n - the number of parameters to be varied (dimension of x).
c
c    a global change of the parameter statements for maxdim & maxsav
c    may be necessary:
c
c    maxdim - the max number of paramters that can be varied.
c             n cannot exceed maxdim.
c    maxsav - the max number of function values that can found (all
c             function values are saved in the tables xsaved, fsaved)
c             maxcl1 and maxcl2 are limited without complaint to
c             maxsav.
c
c     subroutine monitr below is used to print out the progress
c     of the fit.  you may want to modify it.  in addition there
c     are summarizing and error printouts elsewheres in the source.
c
c     8/15/86 - first version - vitaly fiks
c     8/18/86 - small changes - s.p.
c     2/27/86 - small changes - r.b.w.
c     11/27/92 - substituted linpack for nag calls - r.b.w.
c     12/18/92 - changed indexing in quadr - r.b.w.
c
c ----------------------------------------------------------------------
      subroutine minimi(maxcl1,tol1,maxcl2,tol2,a,b,c,scale,
     &                x,fbest,funk,n)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (maxdim=15,maxdm1=maxdim+1,maxsav=200)
      parameter (nlog=0,nin=5,nout=6)
      external funk
      dimension simplx(maxdim,maxdm1),fncval(maxdm1),x(maxdim),
     & scale(maxdim),w1(maxdim),w2(maxdim),
     &  rank(maxdm1),lrank(maxdm1),censim(maxdim),
     & xsaved(maxsav,maxdim),fsaved(maxsav)
c
      if ( n .gt. maxdim ) stop 4444
c
      do 19  i = 1, n
         x(i) = x(i)*scale(i)
 19   continue
      nen = n+1
      maxf = min(maxsav, maxcl1)
      call smplex(simplx,fncval,x,scale,w1,w2,value,n,nen,maxf,tol1
     &            ,rank,ifail,lrank,a,b,c,censim,xsaved,fsaved,ncall,
     &             funk)
      write(nlog,99999) ncall,value,(x(i)/scale(i),i=1,n)
      write(nout,99999) ncall,value,(x(i)/scale(i),i=1,n)
99999 format (/5x,'after ',i3,' function values ',
     &  'final simplex minimum is',f12.5,' at'/5x,(10f12.5)/)
      if (ifail.ne.0) then
        write(nlog,99997) ifail
        write(nout,99997) ifail
99997   format (/5x,'error number:',i3)
      end if
      if (maxcl2.gt.0) then
        maxf = min( maxcl2, maxsav-ncall )
        call quadr(ncall,xsaved,fsaved,maxf,tol2,scale,funk,x,fbest,
     &    n,nen)
      end if
      do 299  i = 1, n
         x(i) = x(i)/scale(i)
 299  continue
      return
      end
c *id* monitr **********************************************************
c ----------------------------------------------------------------------
      subroutine monitr(simplx,scale,fv,fminim,n,nuse,ncall)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (maxdim=15)
      parameter (nlog=0,nin=5,nout=6)
      dimension simplx(maxdim),scale(maxdim)
      if (ncall.eq.1) write(nout,1111)
 1111 format('1')
      write(nlog,9997) ncall,fv,(simplx(i)/scale(i),i=1,n)
      write(nout,9997) ncall,fv,(simplx(i)/scale(i),i=1,n)
 9997 format (' call # ',i3,' gives',f12.5,' at',/5x,(10f12.5))
      if (nuse.ge.0.and.ncall.gt.1) then
        write(nlog,9998) fminim
        write(nout,9998) fminim
 9998   format(' best previous value is',f12.5)
      else if (nuse.lt.0) then
        write(nlog,9999) -nuse
        write(nout,9999) -nuse
 9999   format(' from quadratic fit to ',i3,' rows')
      end if
      return
      end
c *id* functn **********************************************************
c ----------------------------------------------------------------------
      subroutine functn(simp,scale,fv,n,ncall,fsaved,xsaved,funk)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (maxsav=200,maxdim=15)
      external funk
      logical lprt
      dimension simp(n),scale(n),
     1   fsaved(maxsav),xsaved(maxsav,n),xtemp(maxdim)
      ncall=ncall+1
      do 241 i=1,n
        xsaved(ncall,i)=simp(i)
        xtemp(i)=simp(i)/scale(i)
  241 continue
      lprt=.false.
#ifndef ALLOW_TAPENADE
      call funk(xtemp,n,fv,lprt)
#else
      call funk(xtemp,n,fv,ncall)
#endif
      fsaved(ncall)=fv
      return
      end
c *id* smplex **********************************************************
c ----------------------------------------------------------------------
      subroutine smplex(simplx,fncval,x,scale,clexcn,clrefl,value,n,
     &nen,maxcl1,tol1,rank,ifail,lrank,a,b,c,censim,xsaved,fsaved,ncall,
     &funk)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (maxdim=15, maxsav=200)
      external funk
      dimension simplx(maxdim,nen),fncval(nen),
     & scale(n),clexcn(n),clrefl(n)
     &  ,x(n),rank(nen),lrank(nen),censim(n),
     &  fsaved(maxsav),xsaved(maxsav,n)
      ncall=0
      do 1 i=1,n
        simplx(i,1)=x(i)
        do 2 j=2,nen
          simplx(i,j)=simplx(i,1)
    2   continue
        simplx(i,i+1)=simplx(i,i+1)+1
    1 continue
   24 do 12 i=1,nen
        rank(i)=1.e20
   12 continue
      do 3 j=1,nen
        call functn(simplx(1,j),scale,fncval(j),n,ncall,fsaved,xsaved,
     &              funk)
        call monitr(simplx(1,j),scale,fncval(j),rank(nen),n,lrank(nen),
     &              ncall)
        call versrt(fncval(j),lrank,rank,nen,j)
    3 continue
   25 if (ncall.gt.maxcl1) then
         ifail=2
         go to 50
      end if
      call checkt(fncval,nen,chrslt)
      if ((chrslt.lt.tol1).and.(ncall.gt.(nen+1)*nen/2+3*n)) then
         ifail=0
         go to 50
      end if
      call centrd(simplx,censim,n,nen,lrank(1))
      call reflct(simplx,censim,lrank(1),clrefl,n,nen,a)
      call functn(clrefl,scale,fvrefl,n,ncall,fsaved,xsaved,funk)
      call monitr(clrefl,scale,fvrefl,rank(nen),n,lrank(nen),ncall)
      if (fvrefl.lt.rank(nen)) then
         call expand(censim,clrefl,clexcn,n,c)
         call functn(clexcn,scale,fvexpd,n,ncall,fsaved,xsaved,funk)
         call monitr(clexcn,scale,fvexpd,fvrefl,n,ncall-1,ncall)
         if (fvexpd.lt.rank(nen)) then
            call movevt(clexcn,simplx,n,nen,lrank(1),fvexpd,fncval)
            call versrt(fvexpd,lrank,rank,nen,lrank(1))
         end if
         call movevt(clrefl,simplx,n,nen,lrank(1),fvrefl,fncval)
         call versrt(fvrefl,lrank,rank,nen,lrank(1))
      else if (fvrefl.gt.rank(2)) then
         if (fvrefl.lt.rank(1)) then
            call movevt(clrefl,simplx,n,nen,lrank(1),fvrefl,fncval)
            call versrt(fvrefl,lrank,rank,nen,lrank(1))
         end if
         call cntrct(simplx,censim,lrank(1),clexcn,n,nen,b)
         call functn(clexcn,scale,fvcntr,n,ncall,fsaved,xsaved,funk)
         call monitr(clexcn,scale,fvcntr,rank(nen),n,lrank(nen),ncall)
         if (fvcntr.lt.rank(1)) then
            call movevt(clexcn,simplx,n,nen,lrank(1),fvcntr,fncval)
            call versrt(fvcntr,lrank,rank,nen,lrank(1))
         else
            call newset(simplx,fncval,lrank(nen),n,nen)
            go to 24
         end if
      else
         call movevt(clrefl,simplx,n,nen,lrank(1),fvrefl,fncval)
         call versrt(fvrefl,lrank,rank,nen,lrank(1))
      end if
      go to 25
   50 do 60 k=1,n
        x(k)=simplx(k,lrank(nen))
   60 continue
      value=rank(nen)
      return
      end
c *id* versrt **********************************************************
c ----------------------------------------------------------------------
      subroutine versrt(fv,lrank,rank,nen,l)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension rank(nen),lrank(nen)
      lf=l
      i=nen+1
  100 i=i-1
      if (fv.gt.rank(i)) go to 100
      do 120 j=1,i-1
        rank(j)=rank(j+1)
        lrank(j)=lrank(j+1)
  120 continue
      rank(i)=fv
      lrank(i)=lf
      return
      end
c *id* centrd **********************************************************
c ----------------------------------------------------------------------
      subroutine centrd (sim,cen,n,nen,ngr)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (maxdim=15)
      dimension sim(maxdim,nen), cen(n)
      do 9 i=1,n
        cen(i)=0.
    9 continue
      do 10 i=1,n
        do 11 j=1,nen
          if (j.eq.ngr) go to 11
          cen(i)=cen(i)+sim(i,j)
   11   continue
      cen(i)=cen(i)/n
   10 continue
      return
      end
c *id* reflct **********************************************************
c ----------------------------------------------------------------------
      subroutine reflct(simp,cen,ngr,wnewlc,n,nen,alpha)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (maxdim=15)
      dimension simp(maxdim,nen),cen(n),wnewlc(n)
      do 25 j=1,n
        wnewlc(j)=(1+alpha)*cen(j)-alpha*simp(j,ngr)
   25 continue
      return
      end
c *id* expand **********************************************************
c ----------------------------------------------------------------------
      subroutine expand(cen,woldlc,wnewlc,n,gamma)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension cen(n),woldlc(n),wnewlc(n)
      do 26 j=1,n
        wnewlc(j)=gamma*woldlc(j)+(1-gamma)*cen(j)
   26 continue
      return
      end
c *id* cntrct **********************************************************
c ----------------------------------------------------------------------
      subroutine cntrct(simp,cen,ngr,wnewlc,n,nen,beta)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (maxdim=15)
      dimension simp(maxdim,nen),cen(n),wnewlc(n)
      do 27 j=1,n
        wnewlc(j)=beta*simp(j,ngr)+(1-beta)*cen(j)
   27 continue
      return
      end
c *id* movevt **********************************************************
c ----------------------------------------------------------------------
      subroutine movevt(posnew,simplx,n,nen,l,fv,fncval)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (maxdim=15)
      dimension posnew(n),simplx(maxdim,nen),fncval(nen)
      do 100 j=1,n
        simplx(j,l)=posnew(j)
  100 continue
      fncval(l)=fv
      return
      end
c *id* checkt **********************************************************
c ----------------------------------------------------------------------
      subroutine checkt(fncval,nen,sigma)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension fncval(nen)
      fmean=0.
      sigma=0.
      do 16 i=1,nen
        fmean=fmean+fncval(i)
   16 continue
      fmean=fmean/nen
      do 17 i=1,nen
        sigma=sigma+(fncval(i)-fmean)**2
   17 continue
      sigma=sqrt(sigma/nen)
      return
      end
c *id* newset **********************************************************
c ----------------------------------------------------------------------
      subroutine newset(simplx,fncval,nsmlst,n,nen)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      parameter (maxdim=15)
      dimension simplx(maxdim,nen),fncval(n)
      do 35 i=1,n
        do 34 j=1,nen
          simplx(i,j)=simplx(i,nsmlst)+(simplx(i,j)-simplx(i,nsmlst))/2
   34   continue
   35 continue
      return
      end
c *id* quadr ***********************************************************
c ----------------------------------------------------------------------
      subroutine quadr(ncall,xsaved,fsaved,maxcl2,tol2,scale,funk,
     &                  sm,fsmall,n,nen)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      external funk
      parameter (maxdim=15, maxsav=200)
      parameter (nmax=(maxdim+1)*(maxdim+2)/2, maxsin=nmax+15)
      parameter (nlog=0,nin=5,nout=6)
      dimension fsaved(maxsav),xsaved(maxsav,n),scale(n)
      dimension x(maxsin,nmax),s(nmax),nrank(maxsav)
     &,v(nmax,nmax),b(maxsin),c(nmax),d(maxsin)
     &,am(maxdim,maxdim),bm(maxdim),cm(maxdim),sm(maxdim)
     &,ipvt(maxdim),z(maxdim)
      ip=(nen+1)*nen/2
c -------------------------------
c index fsaved in ascending order
c -------------------------------
      ncll=ncall
 7777 call indexx(ncll,fsaved,nrank)
      do 298 i=1,n
        sm(i)=xsaved(nrank(1),i)
  298 continue
      fsmall=fsaved(nrank(1))
      nrow=ip-1
      nfnish=ip
      job=21
      iretry = 0
      ncllsv = ncll
 5678 nrow=nrow+1
      if (nrow.le.ncllsv .and. nrow .le. maxsin) go to 345
      write(nlog,5455)
      write(nout,5455)
 5455 format (/5x,'not enough rows to evaluate the matrix for ',
     &           'the quadratic approximation')
      go to 8120
  345 do 205 i=1,nfnish
        ii = nrank(i)
        x(i,1)=1
        b(i)=fsaved(ii)
        do 206 j=2,nen
          x(i,j)=xsaved(ii,j-1)
  206   continue
        j=nen
        do 207 k=1,n
          do 207 l=k,n
            j=j+1
            x(i,j)=xsaved(ii,k)*xsaved(ii,l)
  207   continue
  205 continue
      if (nrow.ge.ip+5) go to 20
      call dsvdc(x,maxsin,nrow,ip,s,c,x,maxsin,v,nmax,d,job,info)
      do 400 i=1,ip
      if ((s(i)/s(1)).lt.1e-10) go to 401
  400 continue
      go to 20
  401 nfnish=nfnish+1
      go to 5678
   20 do 405 j=1,ip
      c(j)=0.
      do 405 i=1,nrow
        c(j)=x(i,j)*b(i)+c(j)
  405 continue
      do 406 i=1,ip
        d(i)=c(i)/s(i)
  406 continue
      do 407 i=1,ip
      c(i)=0.
      do 407 j=1,ip
        c(i)=d(j)*v(i,j)+c(i)
  407 continue
      lcount=nen
      do 210 i=1,n
        bm(i)=-c(i+1)
        do 210 j=i,n
          lcount=lcount+1
          cc=c(lcount)
          if (i.eq.j) cc=2*cc
          am(i,j)=cc
          am(j,i)=cc
  210 continue
      call dgeco(am,maxdim,n,ipvt,rcond,d)
      call dgesl(am,maxdim,n,ipvt,bm,0)
      do 208 i=1,n
        cm(i)=bm(i)
  208 continue
      if (1.0+rcond.ne.1.0) go to 220
      write(nlog,9995)
      write(nout,9995)
 9995 format (/5x,'the matrix for the derivatives in quadratic
     &             approximation is singular')
      go to 9220
 220  fgrtst=fsaved(nrank(nrow))
      call functn(cm,scale,fqsum,n,ncll,fsaved,xsaved,funk)
      call monitr(cm,scale,fqsum,fqsum,n,-nrow,ncll)
      if (fqsum.lt.fgrtst) go to 521
      write(nlog,8129) fqsum,fgrtst,fsmall
      write(nout,8129) fqsum,fgrtst,fsmall
 8129 format (' new value',f12.5,' is greater than the greatest value',
     & f12.5,' used',/,' best value is',f12.5)
      iretry = iretry + 1
      if ( iretry .lt. 6  .and.  maxcl2 .gt. ncll-ncall) go to 401
      go to 8120
 9220 do 9221  i = 1, ncll
      write(nlog,9222) (xsaved(i,j)/scale(j),j=1,n),fsaved(i)
      write(nout,9222) (xsaved(i,j)/scale(j),j=1,n),fsaved(i)
 9221 continue
 9222 format (4g20.10)
      write(nlog,9223) (s(i),i=1,ip)
      write(nout,9223) (s(i),i=1,ip)
 9223 format (/5x,'singular values are',/(5x,3g20.10))
      write(nlog,9224) (c(i),i=1,ip)
      write(nout,9224) (c(i),i=1,ip)
 9224 format (/5x,'polynomial coef.',/(5x,3g20.10))
      go to 8120
  521 if (abs(fqsum-fsmall).gt.tol2) go to 522
      if ( fsmall .lt. fqsum )  go to 8120
      do 477 i=1,n
        sm(i)=cm(i)
  477 continue
         fsmall=fqsum
 8120 write(nlog,8128) ncll,fsmall,(sm(i)/scale(i),i=1,n)
      write(nout,8128) ncll,fsmall,(sm(i)/scale(i),i=1,n)
 8128 format (/5x,'after ',i3,' function calls final value is',
     & f12.5,' at',/5x,(10f12.5) )
      go to 9881
  522 if (maxcl2.gt.ncll-ncall) go to 7777
      if ( fsmall .lt. fqsum )  go to 8122
      do 577 i=1,n
        sm(i)=cm(i)
  577 continue
         fsmall=fqsum
 8122 write(nlog,8127)
      write(nout,8127)
 8127 format (/5x,'max. number of calls for quad. appr. exceeded')
      go to 8120
 9881 return
      end
c *id* indexx **********************************************************
c indexes an array arrin of length n, ouputting array indx such
c that arrin(indx(j)) is in ascending order for j=1,2,...n.
c arrin is not changed
c see sec 8.3 of numerical recipes by press, et al.
c ----------------------------------------------------------------------
      subroutine indexx(n,arrin,indx)
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension arrin(n),indx(n)
      do 5 j=1,n
        indx(j)=j
    5 continue
      if (n.eq.1) return
      l=n/2+1
      ir=n
   10 continue
      if (l.gt.1) then
        l=l-1
        indxt=indx(l)
        q=arrin(indxt)
      else
        indxt=indx(ir)
        q=arrin(indxt)
        indx(ir)=indx(1)
        ir=ir-1
        if (ir.eq.1) then
          indx(1)=indxt
          return
        end if
      end if
      i=l
      j=l+l
   20 if (j.le.ir) then
        if (j.lt.ir) then
          if (arrin(indx(j)).lt.arrin(indx(j+1))) j=j+1
        end if
        if (q.lt.arrin(indx(j))) then
          indx(i)=indx(j)
          i=j
          j=j+j
        else
          j=ir+1
        end if
        go to 20
      end if
      indx(i)=indxt
      go to 10
      end
