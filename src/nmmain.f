c *id* nmmain **********************************************************
c main subroutine for calculating properties of nuclear or
c neutron matter
c **********************************************************************
      subroutine nmmain(np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
     &                 ,dor,bst,btn,bls,npi,npf, gint, endiff, efree)
      use nmvar
      use nmhncmod
      use nmtbimod
      use nmsubmod
      implicit none
      !implicit real*8 (a-h,o-z)
      !implicit integer*4 (i-n)
      !include "params.f"
      !parameter (ngrid=(20*lgrid+1))
c ----------------------------------------------------------------------
      real*8 :: dor,bst,btn,bls,endiff,efree,gint(6)
      integer*4 :: np,nv,nt,ni,nie,no,ns,lf,lc,ls,lt,ll,lg,le,l3,lk
      integer*4 :: npi,npf,l1,l2
      real*8 xsq(lgrid),xqq(lgrid),rllpp(lgrid),rdls(lgrid),rlss(lgrid)
     &,rlltp(lgrid),wkx(10),wjx(10),ya(6),zif(6,2),esum(14),echeck(14)
     &,vem(14)
      real*8 :: r0,dc,dg,de,ds,dl,dcor,dtor,x,wv2,wvh,wvs,wvs2,wvsb,wfs
      real*8 :: wfs2,wfsb,wj2,wjh,wjs,wjs2,wjsb,wp2,wph,wps,wps2,wpsb
      real*8 :: wk2,wkh,wks,wks2,wksb,wf2,wfh,x1,x2,qx,di1dk1,ditdkt
      real*8 :: wdd,wde,wed,wcc,afe,dli,dlj,dlk,y1,y2,y3,y4,ydd,yde,yee
      real*8 :: ycc,zdd,zde,zee,zcc,yb,ybd,ybc,ybl,z,xdd,xde,xee,xcc,z3
      real*8 :: zk,zf,zj,bvpk,bkpk,bjpk,wv2b,wvhb,wk2b,wf2b,wfhb,wj2b
      real*8 :: wjhb,wp2b,wphb,wv2q,wvhq,wvsq,wvsxq,wkhb,yd,ye,yc,xb
      real*8 :: xbex,xbe,xbc,xi,xicc,zi,bvik,bkik,bjik,bvif,bkif,bjif
      integer*4:: ltd,irl,l,ir,kj,j,i,k,n,mp,m,nx,il,kl,lb,jp,jq,jrl
      real*8 :: bvpf,yt,yr,yl,xcp,zr,zl,fdi,fdk,zkp,zkr,zkl,zkx,bkpf
      real*8 :: bjpf,x3,x4,x6,ytd,ytc,yqd,yqc,elj,yt1,yt2,yt4,xq,xqex
      real*8 :: w3,vemtot,vc1pp,vc1np,vmmpp,vmmnp,vmmnn,esq,tf,wv,wk
      real*8 :: wf,wj,wp,epb,ejf,ev2,evm,ek2,ekm,ej2,ejm,energy
      real*8 :: dlc,dls,dlt,dll,dlg,dle
      dlg=(lg*1.0)
      dle=(le*1.0)
      dlc=(lc*1.0)
      dls=(ls*1.0)
      dlt=(lt*1.0)
      dll=(ll*1.0)
      if (kf.eq.0) kf=(1.5*nm*pi**2*rho)**(1./3.)
      if (rho.eq.0) rho=kf**3/(1.5*nm*pi**2)
      r0=(3/(4*pi*rho))**(1./3.)
      if (nm.eq.1) then
        dt=dor*r0
        dr=dt/(dlt)
      else if (nm.ge.2) then
        dc=dor*r0
        dr=dc/(dlc)
      end if
      ltd=2*lt
      dg=dr*dlg
      de=dr*dle
      dc=dr*dlc
      ds=dr*dls
      dt=dr*dlt
      dl=dr*dll
      dcor=dc/r0
      dtor=dt/r0
c print ================================================================
      if (no.gt.0) then
        write(nlog,1015) kf,rho,r0,dc,ds,dt,dl
        write(nout,1015) kf,rho,r0,dc,ds,dt,dl
 1015   format(/4x,'kf',6x,'rho',5x,'r0',6x,'dc',6x,'ds',6x,'dt',6x,'dl'
     &        ,/7f8.3)
        write(nlog,1016) lc,lt,dcor,dtor,acn,ast,atn,als
        write(nout,1016) lc,lt,dcor,dtor,acn,ast,atn,als
 1016   format(/3x,'lc/lt',3x,'dc/r0',3x,'dt/r0',4x,'ac',6x,'as'
     &         ,6x,'at',6x,'al'/3x,i2,'/',i2,6f8.3)
      end if
c ======================================================================
c -------------------------------------------------------
c find correlation functions and solve fhnc/soc equations
c -------------------------------------------------------
      do 4 irl=1,lgrid
        do jrl=1,8
          gnn(irl,jrl)=0
          v(irl,jrl)=0
        end do
    4 continue
      do 6 irl=1,lgrid
        do jrl=1,8
        f(irl,jrl)=0
        fp(irl,jrl)=0
        fds(irl,jrl)=0
        end do
    6 continue
c ======================================================================
      call nmfts(lc,ls,lt,ll,lf,no,np,nt,nv)
c ======================================================================
      ksav=3*kf**2/5
      kqav=3*kf**4/7
c ======================================================================
      if (temp.ne.0.) call nmhot(lk,ltd,no)
c ======================================================================
      if (ns.eq.1) call nmsps(ltd)
c ======================================================================
      do 13 l=1,4,nm
        if (.not.(l.eq.1)) then
        do 12 ir=1,lt
          f(ir,l)=bst*f(ir,l)
          fp(ir,l)=bst*fp(ir,l)
          fds(ir,l)=bst*fds(ir,l)
   12   continue
        end if
   13 continue
      if (.not.(nv.le.4)) then
      do 14 l=5,6,nm
      do 14 ir=1,lt
        f(ir,l)=btn*f(ir,l)
        fp(ir,l)=btn*fp(ir,l)
        fds(ir,l)=btn*fds(ir,l)
   14 continue
      if (.not.(nv.le.6)) then
      do 16 l=7,8,nm
      do 16 ir=1,lt
        f(ir,l)=bls*f(ir,l)
        fp(ir,l)=bls*fp(ir,l)
        fds(ir,l)=bls*fds(ir,l)
   16 continue
      end if
      end if
c print ================================================================
      if (no.gt.0) then
        write(nlog,1030) dg,de,bst,btn,bls
        write(nout,1030) dg,de,bst,btn,bls
 1030   format(/4x,'dg',6x,'de',6x,'bs',6x,'bt',6x,'bl'/5f8.3)
      end if
c ======================================================================
      call nmhnc(lg,le,l3,ni,nie,no,nt,nv)
c ======================================================================
      do 30 ir=1,ltd
        xsq(ir)=ksav*rs(ir)/3
        xqq(ir)=.16*kqav*rs(ir)**2+2*xsq(ir)
        rllp(ir)=r(ir)*sl(ir)*slp(ir)
        rdls(ir)=rllp(ir)+rs(ir)*(slps(ir)+sl(ir)*sldp(ir))
        rlss(ir)=rllp(ir)+rdls(ir)
        rllpp(ir)=.5*(r(ir)*(slps(ir)+sl(ir)*sldp(ir))-sl(ir)*slp(ir))
        rsdsl(ir)=rs(ir)*sldp(ir)+2*r(ir)*slp(ir)
        rlssx(ir)=gl(ir)*rsdsl(ir)+rs(ir)*slps(ir)
        rlltp(ir)=r(ir)*(sl(ir)*sltp(ir)+3*slp(ir)*sldp(ir))
     &           +2*(slp(ir)**2+sl(ir)*(sldp(ir)-slp(ir)/r(ir)))
   30 continue
c -----------------------
c calculate w0,ws,wf0,wfs
c -----------------------
   40 x=2*rho*pi*dr
      do 50 kj=1,14
        do l1=1,10
          wx(kj,l1)=0
        end do 
   50 continue
      do 51 j=1,10
        wkx(j)=0
        wjx(j)=0
   51 continue
      do 52 kj=1,6
        do l1=1,4
          do l2=1,2
            w3x(kj,l1,l2)=0
          end do
       end do
   52 continue
      wv2=0
      wvh=0
      wvs=0
      wvs2=0
      wvsb=0
      wk2=0
      wkh=0
      wks=0
      wks2=0
      wksb=0
      wf2=0
      wfh=0
      wfs=0
      wfs2=0
      wfsb=0
      wj2=0
      wjh=0
      wjs=0
      wjs2=0
      wjsb=0
      wp2=0
      wph=0
      wps=0
      wps2=0
      wpsb=0
      wx(1,1)=evx
c ---------------
c v6 interactions
c ---------------
      do 150 j=1,6,nm
      do 145 i=1,6,nm
      do 140 k=1,i,nm
      x2=acex(i,j,k)
      if (.not.(x2.eq.0.)) then
      x1=ac(i,j,k)
      qx=2*x
      if (i.eq.k) qx=x
      di1dk1=at(i,1)*at(k,1)
      ditdkt=(at(i,5)+at(i,6))*(at(k,5)+at(k,6))
c -------------------------
c vertex corrections for ws
c -------------------------
      wdd=0
      wde=0
      wed=0
      wcc=0
      do 110 l=1+nm,6,nm
      afe=acex(l,1,l)
      dli=ad(l,i)
      dlj=ad(l,j)
      dlk=ad(l,k)
      y3=.25*(dli+dlk)
      if (.not.(x1.eq.0.)) then
      y1=y3+.25*dlj
      y2=.5*y1+y3/3
      wdd=wdd+y1*bj(l,1)+y2*bj(l,3)
      wde=wde+y3*bj(l,5)
      wed=wed+y1*bj(l,2)+y2*bj(l,4)
      if (i.eq.1.or.j.eq.1.or.k.eq.1) then
        do 95 n=1,4,nm
   95   wde=wde+ak(l,l,n)*aa(n)*(y1+.125*(ad(i,n)+ad(k,n))*at(j,1))
     &   *bj(l,6)/afe
      end if
      end if
      do 105 n=1,4,nm
      y4=.25*(dlj+ad(l,n))
      do 105 mp=1,6,nm
  105 wcc=wcc+.5*(ak(n,i,mp)*ak(j,k,mp)+ak(n,k,mp)*ak(i,j,mp))*aa(mp)
     & *((.5*ad(l,mp)+y4)*bj(l,2)+(.25*ad(l,mp)+.5*y4+y3/3)*bj(l,4))/x2
  110 continue
      ydd=2*(wdd+wde)
      yde=wdd+wde+wed
      yee=2*wed
      ycc=2*wcc
      zdd=(1+wde)**2*exp(2*wdd)-ydd-1
      zde=(1+wde)*exp(wdd+wed)-yde-1
      zee=exp(yee)-yee-1
      zcc=exp(ycc)-ycc-1
c -------------------------------
c vertex corrections for wsb(bpk)
c -------------------------------
      yb=0
      ybd=0
      ybc=0
      do 120 l=7,8,nm
      yb=yb+bj(l,1)
      dlj=ad(l,j)
      ybd=ybd+.25*(ad(l,i)+dlj+ad(l,k))*bj(l,1)
      do 115 n=1,4,nm
      y4=.25*(dlj+ad(l,n))
      do 115 m=1,6,nm
  115 ybc=ybc+.5*(ak(n,i,m)*ak(j,k,m)+ak(n,k,m)*ak(i,j,m))*aa(m)
     & *(.5*ad(l,m)+y4)*bj(l,1)/x2
  120 continue
      ybd=2*ybd
      ybc=2*ybc
      ybl=2*yb
c --------------------------------------
c integrations for w2,wh,ws,ws2,wsb(bpk)
c --------------------------------------
      do 130 ir=1,ltd
      z=f(ir,i)*v(ir,j)*f(ir,k)*qx*rs(ir)
      xdd=gx(ir)
      xde=gy(ir)
      xee=xde**2+gz(ir)
      xcc=gl(ir)**2/nu
      wx(j,1)=wx(j,1)+z*x1
      wx(j,2)=wx(j,2)-z*x2*sls(ir)/nu
      wx(j,3)=wx(j,3)+z*x1*((1+2*xde+xee)*xdd-1)
      wx(j,4)=wx(j,4)-z*x2*(xcc*xdd-sls(ir)/nu)
      wx(j,5)=wx(j,5)+z*x1*(ydd+2*xde*yde+xee*yee)*xdd
      wx(j,6)=wx(j,6)-z*x2*xcc*ycc*xdd
      wx(j,7)=wx(j,7)+z*x1*(zdd+2*xde*zde+xee*zee)*xdd
      wx(j,8)=wx(j,8)-z*x2*xcc*zcc*xdd
      wx(j,9)=wx(j,9)+z*x1*ybd*xdd
      wx(j,10)=wx(j,10)-z*x2*(.5*xcc*ybc+(rlss(ir)*ybl+rs(ir)*slps(ir)
     & *ybc)/(6*nu*xsq(ir)))*xdd
c ----------------------------
c effective 2-body part of tni
c ----------------------------
      z3=f(ir,i)*f(ir,k)*qx*rs(ir)
      w3x(j,1,1)=w3x(j,1,1)
     & +z3*(v3dd(ir,j,1)*(1+ydd+zdd+2*xde*(1+yde+zde)
     & +xee*(1+yee+zee))+2*v3de(ir,j,1)*(1+yde+zde+xde*(1+yee+zee))
     & +v3ee(ir,j,1)*(1+yee+zee))*xdd*x1
      w3x(j,2,1)=w3x(j,2,1)
     & -z3*(v3dd(ir,j,1)*xcc-2*v3cc(ir,j,1)*gl(ir))*(1+ycc+zcc)*xdd*x2
c -----------
c s-wave part
c -----------
      w3x(j,1,2)=w3x(j,1,2)
     & +z3*(v3dd(ir,j,2)*(1+ydd+zdd+2*xde*(1+yde+zde)
     & +xee*(1+yee+zee))+2*v3de(ir,j,2)*(1+yde+zde+xde*(1+yee+zee))
     & +v3ee(ir,j,2)*(1+yee+zee))*xdd*x1
      w3x(j,2,2)=w3x(j,2,2)
     & -z3*(v3dd(ir,j,2)*xcc-2*v3cc(ir,j,2)*gl(ir))*(1+ycc+zcc)*xdd*x2
c ---------------------------
c pair distribution functions
c ---------------------------
      gnn(ir,j)=gnn(ir,j)+f(ir,i)*f(ir,k)*(x1*(1+2*xde+xee+ydd+zdd
     & +2*xde*(yde+zde)+xee*(yee+zee)+ybd)-x2*(xcc+xcc*(ycc+zcc)
     & +.5*xcc*ybc+(rlss(ir)*ybl+rs(ir)*slps(ir)*ybc)/(6*nu*xsq(ir))))
     & *xdd*qx/x
  130 continue
      if (j.eq.1) then
        do 135 ir=1,lt
        xdd=gx(ir)
        xde=gy(ir)
        xee=xde**2+gz(ir)
        xcc=gl(ir)**2/nu
c -------------------
c pb kinetic energies
c -------------------
        zk=-.5*h2m*(f(ir,i)*fds(ir,k)+f(ir,k)*fds(ir,i))*qx*rs(ir)
        wkx(1)=wkx(1)+zk*x1
        wkx(2)=wkx(2)-zk*x2*sls(ir)/nu
        wkx(3)=wkx(3)+zk*x1*((1+2*xde+xee)*xdd-1)
        wkx(4)=wkx(4)-zk*x2*(xcc*xdd-sls(ir)/nu)
        wkx(5)=wkx(5)+zk*x1*(ydd+2*xde*yde+xee*yee)*xdd
        wkx(6)=wkx(6)-zk*x2*xcc*ycc*xdd
        wkx(7)=wkx(7)+zk*x1*(zdd+2*xde*zde+xee*zee)*xdd
        wkx(8)=wkx(8)-zk*x2*xcc*zcc*xdd
        wkx(9)=wkx(9)+zk*x1*ybd*xdd
        wkx(10)=wkx(10)-zk*x2*(.5*xcc*ybc+(rlss(ir)*ybl+rs(ir)*slps(ir)
     &   *ybc)/(6*nu*xsq(ir)))*xdd
        zf=h2m*(f(ir,i)*fp(ir,k)+f(ir,k)*fp(ir,i))*slp(ir)*qx*rs(ir)
        wf2=wf2+zf*x2*sl(ir)/nu
        wfh=wfh+zf*x2*(gl(ir)*xdd-sl(ir))/nu
        wfs=wfs+zf*x2*gl(ir)*xdd*ycc/nu
        wfs2=wfs2+zf*x2*gl(ir)*xdd*zcc/nu
c -------------------
c jf kinetic energies
c -------------------
        zj=h2m*(fp(ir,i)*fp(ir,k)*rs(ir)+6*f(ir,i)*f(ir,k)*ditdkt)*qx
        wjx(1)=wjx(1)+zj*x1
        wjx(2)=wjx(2)-zj*x2*sls(ir)/nu
        wjx(3)=wjx(3)+zj*x1*((1+2*xde+xee)*xdd-1)
        wjx(4)=wjx(4)-zj*x2*(xcc*xdd-sls(ir)/nu)
        wjx(5)=wjx(5)+zj*x1*(ydd+2*xde*yde+xee*yee)*xdd
        wjx(6)=wjx(6)-zj*x2*xcc*ycc*xdd
        wjx(7)=wjx(7)+zj*x1*(zdd+2*xde*zde+xee*zee)*xdd
        wjx(8)=wjx(8)-zj*x2*xcc*zcc*xdd
        wjx(9)=wjx(9)+zj*x1*ybd*xdd
        wjx(10)=wjx(10)-zj*x2*(.5*xcc*ybc+(rlss(ir)*ybl+rs(ir)*slps(ir)
     &   *ybc)/(6*nu*xsq(ir)))*xdd
        wp2=wp2-.5*h2m*x2*(f(ir,i)*f(ir,k)-di1dk1)
     &   *rlss(ir)*qx/nu
        wph=wph-.5*h2m*x2*(f(ir,i)*f(ir,k)*xdd-di1dk1)
     &   *rlssx(ir)*qx/nu
        wps=wps-.5*h2m*x2*(f(ir,i)*f(ir,k)*xdd-di1dk1)
     &   *rlssx(ir)*ycc*qx/nu
        wps2=wps2-.5*h2m*x2*(f(ir,i)*f(ir,k)*xdd-di1dk1)
     &   *rlssx(ir)*zcc*qx/nu
  135   continue
      end if
      end if
  140 continue
  145 continue
      wv2=wv2+wx(j,1)+wx(j,2)
      wvh=wvh+wx(j,3)+wx(j,4)
      wvs=wvs+wx(j,5)+wx(j,6)
      wvs2=wvs2+wx(j,7)+wx(j,8)
      wvsb=wvsb+wx(j,9)+wx(j,10)
  150 continue
      wk2=wk2+wkx(1)+wkx(2)
      wkh=wkh+wkx(3)+wkx(4)
      wks=wks+wkx(5)+wkx(6)
      wks2=wks2+wkx(7)+wkx(8)
      wksb=wksb+wkx(9)+wkx(10)
      wj2=wj2+.5*(wjx(1)+wjx(2)+wkx(1)+wkx(2))
      wjh=wjh+.5*(wjx(3)+wjx(4)+wkx(3)+wkx(4))
      wjs=wjs+.5*(wjx(5)+wjx(6)+wkx(5)+wkx(6))
      wjs2=wjs2+.5*(wjx(7)+wjx(8)+wkx(7)+wkx(8))
      wjsb=wjsb+.5*(wjx(9)+wjx(10)+wkx(9)+wkx(10))
      wph=wph-wp2
      bvpk=wvsb
      bkpk=wksb
      bjpk=wjsb
c -------------------------------------------------
c l.s contributions
c     done with following do loops:
c     j=1,6 i=7,8 k=7,8 gives bcb bsb btb terms
c     j=7,8 i=5,6 k=5,6 gives tbt terms
c     j=7,8 i=7,8 k=1,8 gives bbc bbs bbt bbb terms
c -------------------------------------------------
      wv2b=0
      wvhb=0
      wk2b=0
      wkhb=0
      wf2b=0
      wfhb=0
      wj2b=0
      wjhb=0
      wp2b=0
      wphb=0
      wv2q=0
      wvhq=0
      wvsq=0
      wvsxq=0
      if (.not.(nv.le.6)) then
      nx=nm
      if (nm.eq.4) nx=2
      do 250 j=1,8,nx
      if (.not.(nm.eq.4 .and. j.eq.3)) then
      il=7
      if (j.ge.7) il=5
      do 240 i=il,8,nx
      kl=7
      if (j.ge.7) kl=1
      if (i.le.6) kl=5
      do 240 k=kl,i,nx
      if (.not.(nm.eq.4 .and. k.eq.3)) then
      if (i.le.6) then
        x2=acex(i,j,k)
        x1=ac(i,j,k)
      else
        x2=acl2ex(i,j,k)
        x1=acl2(i,j,k)
      end if
      qx=2*x
      if (i.eq.k) qx=x
c     write(nlog,*) i,j,k,qx,x1,x2
c     write(nout,*) i,j,k,qx,x1,x2
c -------------------------------
c vertex corrections for wsb(bik)
c -------------------------------
      yd=0
      ye=0
      yc=0
      yb=0
      do 210 l=1+nm,6,nm
      dli=ad(l,i)
      dlj=ad(l,j)
      dlk=ad(l,k)
      y3=.25*(dli+dlk)
      y1=y3+.25*dlj
      y2=.5*y1+y3/3
      if (.not.(x1.eq.0.)) then
      yd=yd+y1*bj(l,2)+y2*bj(l,4)
      ye=ye+y3*bj(l,5)
      if (l.le.4) then
        yb=yb+y1*bk(l,2)+bk(l,1)
        if (k.le.6) yb=yb+.5*(y1*(bk(l,3)-bk(l,2))-bk(l,1))
      end if
      end if
      do 205 n=1,4,nm
      y4=.25*(dlj+ad(l,n))
      do 205 mp=1,8,nx
      if (.not.(nm.eq.4 .and. mp.eq.3)) then
      yc=yc+.5*(ak(n,i,mp)*ak(j,k,mp)+ak(n,k,mp)*ak(i,j,mp))*ab(mp)
     & *((.5*ad(l,mp)+y4)*bj(l,2)+(.25*ad(l,mp)+.5*y4+y3/3)*bj(l,4))
      end if
  205 continue
  210 continue
      ydd=2*yd+ye+yb
      yde=2*yd
      ycc=al2(i,j,k)*yc/x2
      if (.not.(i.ge.7)) then
      ydd=ydd+ye-yb
      ycc=-36*yc/x2
      end if
c --------------------------------------
c integrations for w2b,whb,wsb(bik),wf2b
c --------------------------------------
  215 do 230 ir=1,ltd
      z=f(ir,i)*v(ir,j)*f(ir,k)*qx*rs(ir)
      xdd=gx(ir)
      xde=gy(ir)
      xee=gz(ir)
      if (i.ge.7) then
        xb=xsq(ir)
        xbex=rllp(ir)/nu
        xcc=r(ir)*gl(ir)*slp(ir)/nu
        xi=rlss(ir)/nu
        xicc=rlssx(ir)/nu
        if (k.ge.7) then
          xbe=bde(ir,2)/xb
          xbc=bcc(ir,2)
        else
          xbe=(.5*(bde(ir,1)+bde(ir,2))+bde(ir,3))/xb
          xbc=.5*(bcc(ir,1)+bcc(ir,2))+bcc(ir,3)
        end if
c-----------------------------------------------------
c Arya Akmal version uses this instead
c       xbe=bde(ir,1)/xb
c       xbc=bcc(ir,1)
c-----------------------------------------------------
      else
        xb=1
        xbex=sls(ir)/nu
        xcc=gl(ir)**2/nu
        xbe=xde+xde**2+xee
        xbc=0
      end if
      wx(j,1)=wx(j,1)+z*x1*xb
      wx(j,2)=wx(j,2)-z*x2*xbex
      wx(j,3)=wx(j,3)+z*x1*xb*((1+xde+xbe)*xdd-1)
      wx(j,4)=wx(j,4)-z*x2*((xcc+xbc)*xdd-xbex)
      wx(j,5)=wx(j,5)+z*x1*xb*(ydd+xde*yde)*xdd
      wx(j,6)=wx(j,6)-z*x2*xcc*ycc*xdd
      if (j.eq.1) then
c -------------------
c pb kinetic energies
c -------------------
        zk=-.5*h2m*(f(ir,i)*fds(ir,k)+f(ir,k)*fds(ir,i))*qx*rs(ir)
        wkx(1)=wkx(1)+zk*x1*xb
        wkx(2)=wkx(2)-zk*x2*xbex
        wkx(3)=wkx(3)+zk*x1*xb*((1+xde+xbe)*xdd-1)
        wkx(4)=wkx(4)-zk*x2*((xcc+xbc)*xdd-xbex)
        wkx(5)=wkx(5)+zk*x1*xb*(ydd+xde*yde)*xdd
        wkx(6)=wkx(6)-zk*x2*xcc*ycc*xdd
        zf=h2m*(f(ir,i)*fp(ir,k)+f(ir,k)*fp(ir,i))*qx*rs(ir)
        wf2b=wf2b+zf*x2*rllpp(ir)/nu
c -------------------
c jf kinetic energies
c -------------------
        zj=h2m*(fp(ir,i)*fp(ir,k)*rs(ir)+(f(ir,i)*fp(ir,k)
     &   +fp(ir,i)*f(ir,k))*r(ir))*qx
        zi=h2m*f(ir,i)*f(ir,k)*qx
        wjx(1)=wjx(1)+(zj+3*zi)*x1*xb
        wjx(2)=wjx(2)-zj*x2*xbex-zi*x2*xi
        wjx(3)=wjx(3)+(zj+3*zi)*x1*xb*((1+xde+xbe)*xdd-1)
        wjx(4)=wjx(4)-zj*x2*((xcc+xbc)*xdd-xbex)-zi*x2*(xicc*xdd-xi)
        wjx(5)=wjx(5)+(zj+3*zi)*x1*xb*(ydd+xde*yde)*xdd
        wjx(6)=wjx(6)-zj*x2*xcc*ycc*xdd-zi*x2*xicc*ycc*xdd
        wp2b=wp2b-.25*h2m*x2*f(ir,i)*f(ir,k)*rlltp(ir)*rs(ir)*qx/nu
      end if
c ---------------------------
c pair distribution functions
c ---------------------------
      gnn(ir,j)=gnn(ir,j)+f(ir,i)*f(ir,k)*(x1*xb*(1+xde+xbe
     & +ydd+xde*yde)-x2*(xcc+xbc+xcc*ycc))*xdd*qx/x
  230 continue
      end if
  240 continue
      wv2b=wv2b+wx(j,1)+wx(j,2)
      wvhb=wvhb+wx(j,3)+wx(j,4)
      wvsb=wvsb+wx(j,5)+wx(j,6)
      end if
  250 continue
      wv2b=wv2b-wv2
      wvhb=wvhb-wvh
      wvsb=wvsb-wvs
      wk2b=wk2b+wkx(1)+wkx(2)-wk2
      wkhb=wkhb+wkx(3)+wkx(4)-wkh
      wksb=wksb+wkx(5)+wkx(6)-wks
      wj2b=wj2b+.5*(wjx(1)+wjx(2)+wkx(1)+wkx(2))-wj2
      wjhb=wjhb+.5*(wjx(3)+wjx(4)+wkx(3)+wkx(4))-wjh
      wjsb=wjsb+.5*(wjx(5)+wjx(6)+wkx(5)+wkx(6))-wjs
      bvik=wvsb-bvpk
      bkik=wksb-bkpk
      bjik=wjsb-bjpk
c --------
c wsb(bif)
c --------
      bvif=wx(1,5)+wx(1,6)+wx(7,5)+wx(7,6)+wx(8,5)+wx(8,6)
      bkif=wkx(5)+wkx(6)
      bjif=.5*(wjx(5)+wjx(6)+wkx(5)+wkx(6))
      yd=0
      ye=0
      do 252 l=1,6
      ya(l)=0
      zif(l,1)=0
      zif(l,2)=0
  252 continue
      do 260 l=1,6,nm
      do 254 ir=1,lt
  254 ya(l)=ya(l)-aa(l)*f(ir,l)*fds(ir,l)
     & *((1+gy(ir))**2+gz(ir))*gx(ir)*rs(ir)*x
      yd=yd+ya(l)
      do 258 m=1,6,nm
      y1=acex(m,1,l)
      do 256 ir=1,lt
  256 ye=ye+y1*f(ir,m)*fds(ir,l)*gl(ir)**2*gx(ir)*rs(ir)*x/nu
  258 continue
  260 continue
      do 264 l=7,8,nm
      lb=l-6
      if (l.eq.7) then
        ydd=(yd+ye)/3
        ycc=(2*yd+4*(ya(1)+ya(2)))/9
      else if (l.eq.8) then
        ydd=yd+ye
        ycc=(-2*yd-4*ya(1)+8*ya(2)+16*(ya(4)+ya(6))/3)/3
      end if
      do 262 ir=1,lt
      qx=x*rs(ir)
      xdd=((1+gy(ir))**2+gz(ir))*gx(ir)
      xcc=gl(ir)**2*gx(ir)/nu
      zif(1,lb)=zif(1,lb)+.25*v(ir,1)*f(ir,l)**2*rs(ir)*xdd*qx
      zif(2,lb)=zif(2,lb)+.25*v(ir,1)*f(ir,l)**2*rs(ir)*xcc*qx
      zif(3,lb)=zif(3,lb)+f(ir,l)*v(ir,l)*f(ir,1)*rs(ir)*xdd*qx
      zif(4,lb)=zif(4,lb)+f(ir,l)*v(ir,l)*f(ir,1)*rs(ir)*xcc*qx
      zif(5,lb)=zif(5,lb)-.25*h2m*f(ir,l)*fds(ir,l)*rs(ir)*xdd*qx
      zif(6,lb)=zif(6,lb)-.25*h2m*f(ir,l)*fds(ir,l)*rs(ir)*xcc*qx
      gnn(ir,1)=gnn(ir,1)+.25*f(ir,l)**2*rs(ir)*(xdd*ydd-xcc*ycc)
  262 gnn(ir,l)=gnn(ir,l)+f(ir,l)*f(ir,1)*rs(ir)*(xdd*ydd-xcc*ycc)
      wx(1,5)=wx(1,5)+zif(1,lb)*ydd
      wx(1,6)=wx(1,6)-zif(2,lb)*ycc
      wx(l,5)=wx(l,5)+zif(3,lb)*ydd
      wx(l,6)=wx(l,6)-zif(4,lb)*ycc
      wkx(5)=wkx(5)+zif(5,lb)*ydd
      wkx(6)=wkx(6)-zif(6,lb)*ycc
      wjx(5)=wjx(5)+zif(5,lb)*ydd
      wjx(6)=wjx(6)-zif(6,lb)*ycc
  264 continue
      do 266 l=1,3
      ya(l)=0
      zif(l,1)=0
  266 zif(l,2)=0
      do 274 l=7,8,nm
      lb=l-6
      do 270 i=1,6,nm
      do 270 k=1,6,nm
      y1=acex(i,1,k)*(ad(l,i)+.5*(ad(2,i)+ad(2,k)))*aa(lb)/3
      do 268 ir=1,lt
  268 ya(lb)=ya(lb)-y1*f(ir,i)*fp(ir,k)*sl(ir)*slp(ir)
     & *gx(ir)*rs(ir)*x/nu
  270 continue
      do 272 ir=1,lt
      xdd=x*rs(ir)*gx(ir)
      zif(1,lb)=zif(1,lb)+.5*v(ir,1)*f(ir,l)**2*rs(ir)*xdd
      zif(2,lb)=zif(2,lb)+f(ir,l)*v(ir,l)*f(ir,1)*rs(ir)*xdd
      zif(3,lb)=zif(3,lb)-.5*h2m*f(ir,l)*fds(ir,l)*rs(ir)*xdd
      gnn(ir,1)=gnn(ir,1)+.5*f(ir,l)**2*rs(ir)*gx(ir)*ya(lb)
  272 gnn(ir,l)=gnn(ir,l)+f(ir,l)*f(ir,1)*rs(ir)*gx(ir)*ya(lb)
      wx(1,5)=wx(1,5)+zif(1,lb)*ya(lb)
      wx(l,5)=wx(l,5)+zif(2,lb)*ya(lb)
      wkx(5)=wkx(5)+zif(3,lb)*ya(lb)
      wjx(5)=wjx(5)+zif(3,lb)*ya(lb)
  274 continue
      bvif=wx(1,5)+wx(1,6)+wx(7,5)+wx(7,6)+wx(8,5)+wx(8,6)-bvif
      bkif=wkx(5)+wkx(6)-bkif
      bjif=.5*(wjx(5)+wjx(6)+wkx(5)+wkx(6))-bjif
      wvsb=wvsb+bvif
      wksb=wksb+bkif
      wjsb=wjsb+bjif
c --------
c wsb(bpf)
c --------
      bvpf=0
      do 290 j=1,6,nm
        do 288 i=1,6,nm
        do 288 k=1,6,nm
          x2=acex(i,j,k)
          if (.not.(x2.eq.0.)) then
          x1=ac(i,j,k)
          qx=x
          ye=0
          yt=0
          yr=0
          yl=0
          do 278 l=7,8,nm
            y1=6*bj(l,1)/ksav
            ye=ye+(1+ad(l,j))*y1/24
            yt=yt+as(j)*(1+ad(l-6,j))*y1/8
            do 276 n=1,4,nm
            do 276 m=1,6,nm
              y2=(2+ad(l,m)+ad(l,j))*y1/24
              yr=yr+ak(i,j,m)*ak(k,n,m)*aa(m)*y2
              yl=yl+ak(n,j,m)*ak(j,k,m)*aa(m)*y2
  276       continue
  278     continue
          do 280 ir=1,lt
            xdd=qx*gx(ir)
            xcc=xdd*sls(ir)/nu
            xcp=xdd*sl(ir)*slp(ir)/nu
            z=fp(ir,i)*v(ir,j)*fp(ir,k)*rs(ir)
            zr=f(ir,i)*v(ir,j)*fp(ir,k)*rs(ir)
            zl=fp(ir,i)*v(ir,j)*f(ir,k)*rs(ir)
            wx(j,9)=wx(j,9)+z*x1*ye*xdd
            wx(j,10)=wx(j,10)-z*x2*ye*xcc-(zr*yr+zl*yl)*xcp
            gnn(ir,j)=gnn(ir,j)+fp(ir,i)*fp(ir,k)*ye*(x1-x2*sls(ir)/nu)
     &       *gx(ir)-(f(ir,i)*fp(ir,k)*yr+fp(ir,i)*f(ir,k)*yl)*sl(ir)
     &       *slp(ir)*gx(ir)/nu
  280     continue
          if (i.ge.5.and.k.ge.5) then
            do 282 ir=1,lt
              xdd=qx*gx(ir)
              xcc=xdd*sls(ir)/nu
              z=f(ir,i)*v(ir,j)*f(ir,k)
              wx(j,9)=wx(j,9)+z*x1*yt*xdd
              wx(j,10)=wx(j,10)-z*x2*yt*xcc
              gnn(ir,j)=gnn(ir,j)+f(ir,i)*f(ir,k)*yt*(x1-x2*sls(ir)/nu)
     &         *gx(ir)/rs(ir)
  282     continue
          end if
          if (j.eq.1) then
            do 284 ir=1,lt
              xdd=qx*gx(ir)
              xcc=xdd*sls(ir)/nu
              xcp=xdd*sl(ir)*slp(ir)/nu
              fdi=fds(ir,i)-2*fp(ir,i)/r(ir)
              if (i.ge.5) fdi=fdi+6*f(ir,i)/rs(ir)
              fdk=fds(ir,k)-2*fp(ir,k)/r(ir)
              if (k.ge.5) fdk=fdk+6*f(ir,k)/rs(ir)
              zk=fdi*fdk*rs(ir)
              zkp=2*fp(ir,i)*fdk*rs(ir)
              zkr=(f(ir,i)*(sl(ir)*sldp(ir)+slps(ir))+fp(ir,i)*sl(ir)
     &         *slp(ir))*fdk*rs(ir)/nu
              zkl=fp(ir,i)*(fdk+2*fp(ir,k)/r(ir))*rs(ir)
              wkx(9)=wkx(9)+h2m*zk*x1*ye*xdd
              wkx(10)=wkx(10)-h2m*(zk*x2*ye*xcc+zkp*x2*ye*xcp+zkr*yr*xdd
     &         -zkl*yl*xcp)
              wjx(9)=wjx(9)+h2m*zk*x1*ye*xdd
              wjx(10)=wjx(10)-h2m*(zk*x2*ye*xcc+zkp*x2*ye*xcp+zkr*yr*xdd
     &         -zkl*yl*xcp)
              if (i.ge.5.and.k.ge.5) then
                zk=fp(ir,i)*fp(ir,k)
                zkp=f(ir,i)*(fdk+2*fp(ir,k)/r(ir))
                zkr=f(ir,i)*fp(ir,k)
                zkl=fp(ir,i)*f(ir,k)
                zkx=f(ir,i)*f(ir,k)/rs(ir)
                wkx(9)=wkx(9)+h2m*(zk*x1*(yt+6*ye)*xdd+6*zkx*x1*yt*xdd)
                wkx(10)=wkx(10)-h2m*(zk*x2*6*ye*xcc-zkp*x2*yt*xcc
     &           +6*(zkr*yr+zkl*yl)*xcp+6*zkx*x2*yt*xcc)
                wjx(9)=wjx(9)+h2m*(zk*x1*(yt+6*ye)*xdd+6*zkx*x1*yt*xdd)
                wjx(10)=wjx(10)-h2m*(zk*x2*6*ye*xcc-zkp*x2*yt*xcc
     &           +6*(zkr*yr+zkl*yl)*xcp+6*zkx*x2*yt*xcc)
              end if
  284       continue
          end if
          end if
  288   continue
        bvpf=bvpf+wx(j,9)+wx(j,10)
  290 continue
      bvpf=bvpf-bvpk
      bkpf=wkx(9)+wkx(10)-bkpk
      bjpf=.5*(wjx(9)+wjx(10)+wkx(9)+wkx(10))-bjpk
      wvsb=wvsb+bvpf
      wksb=wksb+bkpf
      wjsb=wjsb+bjpf
      if (.not.(nv.le.8)) then
c -----------------------------
c l**2 & (l.s)**2 contributions
c -----------------------------
      do 320 j=9,14,nm
      jp=j-8
      jq=j-8
      if (j.ge.13) jq=j-12
      do 316 i=1,8,nx
      if (.not.(nm.eq.4 .and. i.eq.3)) then
      if (j.ge.13.and.i.le.4) jp=j-10
      if (j.ge.13.and.i.ge.5) jp=j-8
      do 315 k=1,i,nx
      if (.not.(nm.eq.4 .and. k.eq.3)) then
      x2=acl2ex(i,j,k)
      if (.not.(x2.eq.0.)) then
      x1=acl2(i,j,k)
      x3=ac(i,j,k)
      x4=acex(i,j,k)
      x6=acex(i,jp,k)
      qx=2*x
      if (i.eq.k) qx=x
c -----------------------------
c passive f6 vertex corrections
c -----------------------------
      yd=0
      yc=0
      ytd=0
      ytc=0
      yqd=0
      yqc=0
      if (.not.(i.ge.7)) then
      do 303 l=1+nm,6,nm
      dli=ad(l,i)
      dlj=ad(l,jp)
      dlk=ad(l,k)
      elj=0
      if (j.ge.13) elj=ae(l,j-12)
      y3=.25*(dli+dlk)
      if (x1.ne.0.) then
        y1=y3+.25*(dlj+elj*acl2(i,j-4,k)/x1)
        y2=.5*y1+y3/3
        yd=yd+y1*bj(l,1)+y2*bj(l,3)+y3*bj(l,5)
      end if
      if (x3.ne.0.) then
        yt1=y3+.25*(dlj+elj*ac(i,j-4,k)/x3)
        yt2=.5*yt1+y3/3
        ytd=ytd+yt1*bj(l,1)+yt2*bj(l,3)+y3*bj(l,5)
      end if
      do 300 n=1,4,nm
      y4=.25*(dlj+elj*acl2ex(i,j-4,k)/x2+ad(l,n))
      do 300 m=1,6,nm
  300 yc=yc+.5*(ak(n,i,m)*ak(jp,k,m)+ak(n,k,m)*ak(i,jp,m))*aa(m)
     & *((.5*ad(l,m)+y4)*bj(l,2)+(.25*ad(l,m)+.5*y4+y3/3)*bj(l,4))/x6
      if (x4.ne.0.) then
        do 302 n=1,4,nm
        yt4=.25*(dlj+elj*acex(i,j-4,k)/x4+ad(l,n))
        do 302 m=1,6,nm
  302   ytc=ytc+.5*(ak(n,i,m)*ak(jp,k,m)+ak(n,k,m)*ak(i,jp,m))*aa(m)
     &   *((.5*ad(l,m)+yt4)*bj(l,2)+(.25*ad(l,m)+.5*yt4+y3/3)*bj(l,4))
     &   /x6
      end if
  303 continue
      yd=2*yd
      yc=2*yc
      ytd=2*ytd
      ytc=2*ytc
c ---------------------------------
c interacting f6 vertex corrections
c ---------------------------------
      do 308 l=1,6,nm
      dli=ad(l,i)
      dlj=ad(l,jq)
      dlk=ad(l,k)
      if (x1.ne.0.) then
        y1=(1+.25*(dli+dlj+dlk))*ac(i,jq,k)
        y2=(1+.25*(dli+dlk))*ac(i,jq,k)
        if (j.ge.13) then
          y1=.5*y1+(.5+.125*(dli+ad(l,jq+2)+dlk))*ac(i,jq+2,k)/3
          y2=.5*y2+(.5+.125*(dli+dlk))*ac(i,jq+2,k)/3
        end if
        yqd=yqd+y1*bq(l,1)+y2*bq(l,2)
      end if
      do 304 n=1,4,nm
      do 304 m=1,6,nm
  304 yqc=yqc+.5*(ak(n,i,m)*ak(jq,k,m)+ak(n,k,m)*ak(i,jq,m))*aa(m)
     & *(1+.5*ad(l,m)+.25*(dlj+ad(l,n)))*bq(l,1)
      if (j.ge.13) then
        yqc=.5*yqc
        do 305 n=1,4,nm
        do 305 m=1,6,nm
  305   yqc=yqc+.5*(ak(n,i,m)*ak(jq+2,k,m)+ak(n,k,m)*ak(i,jq+2,m))*aa(m)
     &   *(.5+.25*ad(l,m)+.125*(ad(l,jq+2)+ad(l,n)))*bq(l,1)/3
      end if
  308 continue
      end if
c ------------------------
c integration for w2q, wsq
c ------------------------
      do 310 ir=1,ltd
      z=f(ir,i)*v(ir,j)*f(ir,k)*qx*rs(ir)
      if (k.le.6) then
        xq=xsq(ir)
        xqex=rllp(ir)/nu
      else
        xq=xqq(ir)
        xqex=rdls(ir)/nu
      end if
      wx(j,1)=wx(j,1)+z*(x1*xq+x3)
      wx(j,2)=wx(j,2)-z*(x2*xqex+x4*sls(ir)/nu)
      wx(j,5)=wx(j,5)+z*(x1*xq*yd+x3*ytd)
      wx(j,6)=wx(j,6)-z*(x2*xqex*yc+x4*ytc*sls(ir)/nu)
      wx(j,9)=wx(j,9)+.5*z*rs(ir)*yqd
      wx(j,10)=wx(j,10)-.5*z*rs(ir)*yqc*sls(ir)/nu
c ---------------------------
c pair distribution functions
c ---------------------------
      gnn(ir,j)=gnn(ir,j)+f(ir,i)*f(ir,k)*(x1*xq*(1+yd)-x2*xqex*(1+yc)
     & +x3*(1+ytd)-x4*(1+ytc)*sls(ir)/nu+.5*rs(ir)*(yqd-yqc*sls(ir)/nu))
     & *qx/x
  310 continue
      end if
      end if
  315 continue
      end if
  316 continue
      wv2q=wv2q+wx(j,1)+wx(j,2)
      wvsq=wvsq+wx(j,5)+wx(j,6)
      wvsxq=wvsxq+wx(j,9)+wx(j,10)
  320 continue
      end if
      end if
      do 410 j=1,10
      wjx(j)=.5*(wkx(j)+wjx(j))
  410 continue
c print ================================================================
      if (.not.(no.eq.0)) then
      write(nlog,1100) wv2,wk2,wf2,wj2,wp2,wv2b,wk2b,wf2b,wj2b,wp2b
      write(nout,1100) wv2,wk2,wf2,wj2,wp2,wv2b,wk2b,wf2b,wj2b,wp2b
 1100 format(/4x,'wv2',5x,'wk2',5x,'wf2',5x,'wj2',5x,'wp2'
     &,5x,'wv2b',4x,'wk2b',4x,'wf2b',4x,'wj2b',4x,'wp2b'/10f8.3)
      write(nlog,1105) wvh,wkh,wfh,wjh,wph,wvhb,wkhb,wfhb,wjhb,wphb
      write(nout,1105) wvh,wkh,wfh,wjh,wph,wvhb,wkhb,wfhb,wjhb,wphb
 1105 format(/4x,'wvh',5x,'wkh',5x,'wfh',5x,'wjh',5x,'wph'
     &,5x,'wvhb',4x,'wkhb',4x,'wfhb',4x,'wjhb',4x,'wphb'/10f8.3)
      write(nlog,1110) wvs,wks,wfs,wjs,wps,wvsb,wksb,wfsb,wjsb,wpsb
      write(nout,1110) wvs,wks,wfs,wjs,wps,wvsb,wksb,wfsb,wjsb,wpsb
 1110 format(/4x,'wvs',5x,'wks',5x,'wfs',5x,'wjs',5x,'wps'
     &,5x,'wvsb',4x,'wksb',4x,'wfsb',4x,'wjsb',4x,'wpsb'/10f8.3)
      write(nlog,1115) wvs2,wks2,wfs2,wjs2,wps2,wv2q,wvhq,wvsq,wvsxq
      write(nout,1115) wvs2,wks2,wfs2,wjs2,wps2,wv2q,wvhq,wvsq,wvsxq
 1115 format(/4x,'wvs2',4x,'wks2',4x,'wfs2',4x,'wjs2',4x,'wps2'
     &,4x,'wv2q',4x,'wvhq',4x,'wvsq',4x,'wvsxq'/9f8.3)
      if (nv.ge.7.and.no.ge.2) then
        write(nout,1120) bvik,bkik,bjik,bvif,bkif,bjif
 1120   format(/4x,'bvik',4x,'bkik',4x,'bjik',4x,'bvif',4x,'bkif',4x
     &  ,'bjif*'/6f8.3)
        write(nout,1125) bvpk,bkpk,bjpk,bvpf,bkpf,bjpf
 1125   format(/4x,'bvpk',4x,'bkpk',4x,'bjpk',4x,'bvpf',4x,'bkpf',4x
     &  ,'bjpf*'/6f8.3)
      end if
      write(nout,1130)
 1130 format(/4x,'w2:d,e')
      write(nout,1133)
 1133 format(4x,'c',7x,'t',7x,'s',7x,'st',6x,'tn',6x,'tnt'
     &,5x,'b',7x,'bt',6x,'ke',6x,'jf')
      write(nout,1140) ((wx(i,j),i=1,8),wkx(j),wjx(j),j=1,2)
 1140 format(10f8.3)
      if (nv.gt.8) then
        write(nout,1130)
        write(nout,1144)
 1144   format(4x,'q',7x,'qt',6x,'qs',6x,'qst',5x,'bb',6x,'bbt')
        write(nout,1150) ((wx(i,j),i=9,14),j=1,2)
 1150   format(6f8.3)
      end if
      write(nout,1160)
 1160 format(/4x,'wh:d,e')
      write(nout,1133)
      write(nout,1140) ((wx(i,j),i=1,8),wkx(j),wjx(j),j=3,4)
      write(nout,1170)
 1170 format(/4x,'ws:d,e')
      write(nout,1133)
      write(nout,1140) ((wx(i,j),i=1,8),wkx(j),wjx(j),j=5,6)
      if (nv.ge.9) then
        write(nout,1170)
        write(nout,1144)
        write(nout,1150) ((wx(i,j),i=9,14),j=5,6)
      end if
      write(nout,1175)
 1175 format(/4x,'ws2:d,e')
      write(nout,1133)
      write(nout,1140) ((wx(i,j),i=1,8),wkx(j),wjx(j),j=7,8)
      if (nv.ge.7) then
        write(nout,1180)
 1180   format(/4x,'wsb:d,e')
        write(nout,1133)
        write(nout,1140) ((wx(i,j),i=1,8),wkx(j),wjx(j),j=9,10)
      end if
      if (nv.ge.9) then
        write(nout,1190)
 1190   format(/4x,'wsx:d,e')
        write(nout,1144)
        write(nout,1150) ((wx(i,j),i=9,14),j=9,10)
      end if
      end if
c ======================================================================
      ev6=wv2+wvh+wvs+wvs2
      evb=wv2b+wvhb+wvsb
      evq=wv2q+wvhq+wvsq+wvsxq
      ek6=wk2+wkh+wks+wks2
      ekb=wk2b+wkhb+wksb
      ef6=wf2+wfh+wfs+wfs2
      ej6=wj2+wjh+wjs+wjs2
      ejb=wj2b+wjhb+wjsb
      ep6=wp2+wph+wps+wps2
c -------------------------
c calculate wc,wcs,wfc,wfcs
c ================================
      call nmchain(nv,nt,no,lg,l3)
c ================================
c calculate u and uf and 3-body part of tni
c =============================================
      call nmtbi(lt,lg,le,l3,nie,no,nt)
c =============================================
      if (nt.eq.4) w3va=-tnia*rho**2*exp(-tnic*rho)*(s-1)
      w3=w3v0+w3v1+w3va+w3vc
c ---------------------
c electromagnetic terms
c ---------------------
      vemtot=0
      if (np.eq.9) then
        vc1pp=0
        vc1np=0
        vmmpp=0
        vmmnp=0
        vmmnn=0
        esq=197.327053/137.03599
        do 750 ir=1,ltd
          call empot(1,1.,1.,r(ir),vem)
          vc1pp=vc1pp+vem(1)*(gnn(ir,1)+gnn(ir,2)/3.)*rs(ir)-esq*r(ir)
          vc1np=vc1np+vem(5)*(gnn(ir,1)-gnn(ir,2)/3.)*rs(ir)
          vmmpp=vmmpp+(vem(6)*(gnn(ir,3)+gnn(ir,4)/3.)
     &                +vem(9)*(gnn(ir,5)+gnn(ir,6)/3.)
     &               +vem(12)*(gnn(ir,7)+gnn(ir,8)/3.))*rs(ir)
          vmmnp=vmmnp+(vem(8)*(gnn(ir,3)-gnn(ir,4)/3.)
     &               +vem(11)*(gnn(ir,5)-gnn(ir,6)/3.)
     &               +vem(14)*(gnn(ir,7)-gnn(ir,8)/3.))*rs(ir)
          vmmnn=vmmnn+(vem(7)*(gnn(ir,3)+gnn(ir,4)/3.)
     &               +vem(10)*(gnn(ir,5)+gnn(ir,6)/3.))*rs(ir)
  750   continue
        if (nm.eq.1) then
          vc1pp=.25*x*vc1pp
          vc1np=.5*x*vc1np
          vmmpp=.25*x*vmmpp
          vmmnp=.5*x*vmmnp
          vmmnn=.25*x*vmmnn
        else if (nm.eq.2) then
          vc1pp=0
          vc1np=0
          vmmpp=0
          vmmnp=0
          vmmnn=x*vmmnn
        end if
        vemtot=vc1pp+vc1np+vmmpp+vmmnp+vmmnn
        if (no.ge.1) then
          write(nlog,1395) vemtot,vc1pp,vc1np,vmmpp,vmmnp,vmmnn
          write(nout,1395) vemtot,vc1pp,vc1np,vmmpp,vmmnp,vmmnn
 1395     format(/4x,'vem',5x,'vc1pp',3x,'vc1np',3x,'vmmpp',3x,'vmmnp'
     &           ,3x,'vmmnn'/6f8.3)
        end if
      end if
c --------
c checksum
c --------
      do 820 l=1,6,nm
      echeck(l)=0
      gint(l)=0
      do 810 ir=1,ltd
      echeck(l)=echeck(l)+v(ir,l)*gnn(ir,l)*rs(ir)*x
      gnn(ir,l)=gnn(ir,l)/aa(l)
      gint(l)=gint(l)+2*(gnn(ir,l)-at(l,1)+af(l)*sls(ir)/nu)*rs(ir)*x
  810 continue
      gint(l)=gint(l)-af(l)
      esum(l)=0
      do 815 j=1,10
      esum(l)=esum(l)+wx(l,j)+wcx(l,j)+wcdx(l,j)+wcmx(l,j)+wcrx(l,j)
  815 continue
  820 continue
      if (nv.gt.6) then
        do 840 l=7,14,nm
        echeck(l)=0
        do 830 ir=1,ltd
        echeck(l)=echeck(l)+v(ir,l)*gnn(ir,l)*rs(ir)*x
  830   continue
        esum(l)=0
        do 835 j=1,10
        esum(l)=esum(l)+wx(l,j)
        if (l.le.8) esum(l)=esum(l)+wcx(l,j)
  835   continue
  840   continue
      end if
c -----------------
c energy summations
c -----------------
      tf=.5*h2m*ksav
      wv=ev6+evb+evq+vemtot
      wk=ek6+ekb
      wf=ef6+wf2b
      wj=ej6+ejb
      wp=ep6+wp2b
      epb=tf+wk+wf+u+uf+wv+w3
      ejf=tf+wj+wp+up+wv+w3
      ev2=wv2+wv2b+wv2q
      evm=wv-ev2
      ek2=wk2+wk2b+wf2+wf2b
      ekm=wk+wf+u+uf-ek2
      ej2=wj2+wj2b+wp2+wp2b
      ejm=wj+wp+up-ej2
      endiff=epb-ejf
      energy=.5*(epb+ejf)
      efree=energy-temp*entrpy
c print ================================================================
      if (no.eq.0) then
      write(nlog,1410) gint
      write(nout,1410) gint
      write(nlog,1430) energy,endiff,epb,ejf
      write(nout,1430) energy,endiff,epb,ejf
      else
      write(nlog,1400) tf,wv,wk,wf,wj,wp,w3,u,uf,up
      write(nout,1400) tf,wv,wk,wf,wj,wp,w3,u,uf,up
 1400 format(/4x,'tf',6x,'wv',6x,'wk',6x,'wf',6x,'wj',6x,'wp',6x,'w3'
     $,6x,'u',7x,'uf',6x,'up'/10f8.3)
      if (no.ge.2) then
        write(nout,1402)
 1402   format(/4x,'esum')
        write(nout,1404)
 1404   format(4x,'c',7x,'t',7x,'s',7x,'st',6x,'tn',6x,'tnt'
     &         ,5x,'b',7x,'bt')
        write(nout,1140) (esum(l),l=1,8)
        write(nout,1402)
        write(nout,1406)
 1406   format(4x,'q',7x,'qt',6x,'qs',6x,'qst',5x,'bb',6x,'bbt')
        write(nout,1140) (esum(l),l=9,14)
        write(nout,1407)
 1407   format(/4x,'echeck')
        write(nout,1404)
        write(nout,1140) (echeck(l),l=1,8)
        write(nout,1407)
        write(nout,1406)
        write(nout,1140) (echeck(l),l=9,14)
      end if
      write(nlog,1410) gint
      write(nout,1410) gint
 1410 format(/4x,'gint'/6f8.3)
      write(nlog,1420) ev6,ek6,ef6,ej6,ep6,evb,ekb,ejb,evq
      write(nout,1420) ev6,ek6,ef6,ej6,ep6,evb,ekb,ejb,evq
 1420 format(/4x,'ev6',5x,'ek6',5x,'ef6',5x,'ej6',5x,'ep6',5x,'evb',5x
     &,'ekb',5x,'ejb',5x,'evq'/9f8.3)
      write(nlog,1425) ev2,evm,ek2,ekm,ej2,ejm
      write(nout,1425) ev2,evm,ek2,ekm,ej2,ejm
 1425 format(/4x,'v(2b)',3x,'v(mb)',3x,'k(2b)',3x,'k(mb)',3x,'j(2b)'
     &,3x,'j(mb)'/6f8.3)
      write(nlog,1430) energy,endiff,epb,ejf
      write(nout,1430) energy,endiff,epb,ejf
 1430 format(/2x,'energy  endiff    epb     ejf',/4f8.3)
      if (temp.ne.0.) then
        write(nlog,1435) efree
        write(nout,1435) efree
 1435   format(/4x,'free energy'/f8.3)
      end if
      end if
c ======================================================================
  850 if (npi.gt.0) call nmpion(ltd,no,np,npi,npf)
      return
      end
