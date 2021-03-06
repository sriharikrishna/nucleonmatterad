c *id* nmhot ***********************************************************
c subroutine for finite temperature
c single-particle spectrum e(k) and occupation number n(k) parameterized
c by temp and chmpot: e(k)=(hbar*k)**2/(2*mstar)
c                     n(k)=sum[ 1/{1+exp[(e(k)-chmpot)/(k*temp)]} ]
c integrals over k from 0 to lkf*kf using simpson's rule
c chmpot found by iterative solution of sum[n(k)]=rho*volume
c entrpy, ksav, kqav, and slater functions calculated
c **********************************************************************
      subroutine nmhot(lkf,ld,no)
      use nmvar
      implicit none
      !implicit real*8 (a-h,o-z)
      !implicit integer*4 (i-n)
      !include "params.f"
      !parameter (ngrid=(20*lgrid+1))

      !real*8 kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,h2m,h2mcs,pi,s
      !common /consts/ kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,
      !&       h2m,h2mcs,pi,s
      !real*8 r(lgrid),ri(lgrid),rs(lgrid),sl(lgrid),sls(lgrid),
      !&       slp(lgrid),slps(lgrid),sldp(lgrid),sltp(lgrid),
      !&       rllp(lgrid),rlssx(lgrid),rsdsl(lgrid)
      !common /rslate/ r,ri,rs,sl,sls,slp,slps,sldp,sltp,rllp,rlssx,rsdsl
      !real*8 temp,mstar,chmpot,entrpy,ksav,kqav
      !common /hotted/ temp,mstar,chmpot,entrpy,ksav,kqav
      !real*8 rx(ngrid),slx(ngrid),slpx(ngrid),sldpx(ngrid),sltpx(ngrid)
      !common /hotfun/ rx,slx,slpx,sldpx,sltpx
      integer*4 :: lkf,ld,no,lk,ik,ii,il,im,in,ir,havetobreak
      real*8 ek(401),k(401),ks(401),nk(401),wk(401),kinc
      real*8 :: sml,htms,fact,tkfc,drho,x,drhos,chmpts,y
      real*8 :: zlx,zlpx,zldpx,zltpx,z
      real*8 :: dlk
c -------------------------------,
c set k, k**2, e(k), weights w(k)
c -------------------------------
      sml=1.e-7
      htms=.5*h2m/mstar
      chmpot=htms*kf**2
      lk=400
      dlk = lk
      kinc=lkf*kf/dlk
      k(1)=0.
      do 5 ik=2,lk+1
        k(ik)=k(ik-1)+kinc
        ks(ik)=k(ik)**2
        ek(ik)=htms*ks(ik)
    5 continue
      wk(1)=1./3.
      wk(lk+1)=wk(1)
      do 6 ik=2,lk,2
        wk(ik)=4./3.
    6 continue
      do 7 ik=3,lk-1,2
        wk(ik)=2./3.
    7 continue
      fact=s*kinc/(2.*pi**2)
      tkfc=3.*kinc/kf**3
c -----------------------
c find chemical potential
c -----------------------
      havetobreak = 0
      do 20 ii=1,20
        if(havetobreak.eq.0) then
        drho=0.
        do 10 ik=1,lk+1
          x=(ek(ik)-chmpot)/temp
          if (x.gt.100.) x=100.
          if (x.lt.-100.) x=-100.
          drho=drho+wk(ik)*ks(ik)/(1.+exp(x))
   10   continue
        drho=fact*drho-rho
        if (.not.(ii.gt.1)) then
        drhos=drho
        chmpts=chmpot
        chmpot=.99*chmpot
        else
   15   if (abs(drhos-drho).lt.sml) then
          havetobreak = 1
        else
        x=(drhos*chmpot-drho*chmpts)/(drhos-drho)
        drhos=drho
        chmpts=chmpot
        chmpot=x
        end if
        end if
        end if
   20 continue
   25 continue
c ---------------------------------
c calculate entropy, <k**2>, <k**4>
c ---------------------------------
      do 30 ik=1,lk+1
        x=(ek(ik)-chmpot)/temp
        if (x.gt.100.) x=100.
        if (x.lt.-100.) x=-100.
        nk(ik)=1./(1.+exp(x))
   30 continue
      entrpy=0.
      ksav=0.
      kqav=0.
      fact=fact/rho
      do 35 ik=1,lk+1
        x=1.-nk(ik)
        y=nk(ik)
        if (x.gt.sml.and.y.gt.sml)
     &    entrpy=entrpy+wk(ik)*ks(ik)*(x*log(x)+y*log(y))
        ksav=ksav+wk(ik)*ks(ik)**2*y
        kqav=kqav+wk(ik)*ks(ik)**3*y
   35 continue
      entrpy=-fact*entrpy
      ksav=fact*ksav
      kqav=fact*kqav
c print ================================================================
      if (.not.(no.eq.0)) then
      write(nlog,1000) mstar,chmpot,lk,lkf,ksav,kqav,entrpy
      write(nout,1000) mstar,chmpot,lk,lkf,ksav,kqav,entrpy
 1000 format(/4x,'finite temperature spectrum has m*/m =',f5.2,
     & ', chemical potential =',f8.3,
     & /4x,'k integration uses',i4,' points from 0 to',i2,' kf'
     & /2x,'<k**2>',2x,'<k**4>',2x,'entropy'/3f8.3)
      if (no.eq.3) then
        write(nout,1010)
 1010   format(/3x,4('k',8x,'n(k)',5x))
        do 37 ik=1,lk/4,5
          il=ik+lk/4
          im=il+lk/4
          in=im+lk/4
          write(nout,1015) k(ik),nk(ik),k(il),nk(il),k(im),nk(im)
     &                    ,k(in),nk(in)
   37   continue
 1015   format(1x,1p,8e9.2)
      end if
      end if
c ======================================================================
c -------------------------------
c calculate slater function, etc.
c -------------------------------
      do 45 ir=1,ld
        zlx=0.
        zlpx=0.
        zldpx=0.
        zltpx=0.
        do 40 ik=1,lk+1
          x=k(ik)*rx(ir)
          y=sin(x)
          z=cos(x)
          zlx=zlx+wk(ik)*nk(ik)*k(ik)*y
          zlpx=zlpx+wk(ik)*nk(ik)*ks(ik)*z
          zldpx=zldpx+wk(ik)*nk(ik)*ks(ik)*k(ik)*y
          zltpx=zltpx+wk(ik)*nk(ik)*ks(ik)**2*z
   40   continue
        slx(ir)=tkfc*zlx/rx(ir)
        slpx(ir)=tkfc*(zlpx-zlx/rx(ir))/rx(ir)
        sldpx(ir)=-tkfc*(zldpx+2.*(zlpx-zlx/rx(ir))/rx(ir))/rx(ir)
        sltpx(ir)=-tkfc*(zltpx-3.*(zldpx+2.*(zlpx-zlx/rx(ir))/rx(ir))
     &                  /rx(ir))/rx(ir)
   45 continue
      return
      end
c *id* nmsps ***********************************************************
c subroutine for single particle spectrum
c **********************************************************************
      subroutine nmsps(ltd)
      use nmvar
      implicit none
      !implicit real*8 (a-h,o-z)
      !implicit integer*4 (i-n)
      !include "params.f"
      !parameter (ngrid=(20*lgrid+1))
      !parameter (nlog=0,nin=5,nout=6)
      !real*8 kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,h2m,h2mcs,pi,s
      !common /consts/ kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,
      !&       h2m,h2mcs,pi,s
      !real*8 r(lgrid),ri(lgrid),rs(lgrid),sl(lgrid),sls(lgrid),
      ! &       slp(lgrid),slps(lgrid),sldp(lgrid),sltp(lgrid),
      ! &       rllp(lgrid),rlssx(lgrid),rsdsl(lgrid)
      !common /rslate/ r,ri,rs,sl,sls,slp,slps,sldp,sltp,rllp,rlssx,rsdsl
      !real*8 temp,mstar,chmpot,entrpy,ksav,kqav
      !common /hotted/ temp,mstar,chmpot,entrpy,ksav,kqav
      !real*8 xph,yph
      !common /parhol/ xph,yph
      integer*4 :: ltd,ir
      real*8 :: rk,ry,dlph,dlphp,dlphdp,dlphtp,dlkf,dlkfp,dlkfdp,dlkftp
c
      ksav=3*kf**2/5+xph*(yph**2-kf**2)
      kqav=3*kf**4/7+xph*(yph**4-kf**4)
c    -------------------------------
c    calculate slater function, etc.
c    -------------------------------
      do 10 ir=1,ltd
      rk=r(ir)*kf
      ry=r(ir)*yph
      dlph=sin(ry)/ry
      dlphp=(cos(ry)-dlph)/r(ir)
      dlphdp=-yph**2*dlph-2*dlphp/r(ir)
      dlphtp=(2/rs(ir)-yph**2)*dlphp-2*dlphdp/r(ir)
      dlkf=sin(rk)/rk
      dlkfp=(cos(rk)-dlkf)/r(ir)
      dlkfdp=-kf**2*dlkf-2*dlkfp/r(ir)
      dlkftp=(2/rs(ir)-kf**2)*dlkfp-2*dlkfdp/r(ir)
      sl(ir)=sl(ir)+xph*(dlph-dlkf)
      sls(ir)=sl(ir)**2
      slp(ir)=slp(ir)+xph*(dlphp-dlkfp)
      slps(ir)=slp(ir)**2
      sldp(ir)=sldp(ir)+xph*(dlphdp-dlkfdp)
      sltp(ir)=sltp(ir)+xph*(dlphtp-dlkftp)
   10 continue
      return
      end

c *id* nmpion **********************************************************
c pion density subroutine
c pion excess npi(q) in %/nucleon
c calculated for q=grid( qmin, qmax, (qmax-qmin)/(npi-1) )
c abs(npf)=1(2) monopole (dipole) form factor
c sign(npf)=+(-) configuration (momentum) wave function
c eav=energy denominator
c fsof=(fpind/fpinn)**2
c plm=lambda (fm**-1) in form factor
c **********************************************************************
      subroutine nmpion(lp,no,np,npi,npf)
      use nmvar
      implicit none
      !implicit real*8 (a-h,o-z)
      !implicit integer*4 (i-n)
      !include "params.f"
      !parameter (ngrid=(20*lgrid+1))
      !parameter (nlog=0,nin=5,nout=6)
      !real*8 kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,h2m,h2mcs,pi,s
      !common /consts/ kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,
      !&       h2m,h2mcs,pi,s
      !real*8 r(lgrid),ri(lgrid),rs(lgrid),sl(lgrid),sls(lgrid),
      ! &       slp(lgrid),slps(lgrid),sldp(lgrid),sltp(lgrid),
      !&       rllp(lgrid),rlssx(lgrid),rsdsl(lgrid)
      !common /rslate/ r,ri,rs,sl,sls,slp,slps,sldp,sltp,rllp,rlssx,rsdsl
      !real*8 aa(8),ab(8),ad(8,8),ae(6,2),af(8),ak(8,8,8),al(6,6,6),
      !&       as(6),at(8,8),ax(6,6,6)
      !common /amatrx/ aa,ab,ad,ae,af,ak,al,as,at,ax
      !real*8 gca(lgrid,6),gcb(lgrid,6),gdd(lgrid,6),gde(lgrid,6),
      !&       gee(lgrid,6),gl(lgrid),gx(lgrid),gy(lgrid),gz(lgrid),
      !&       gnn(lgrid,14)
      !common /gchain/ gca,gcb,gdd,gde,gee,gl,gx,gy,gz,gnn
      !real*8 tpi(lgrid),ypi(lgrid),tpi2(lgrid),
       !&       xt0(lgrid),xt1(lgrid),xt2(lgrid),xt3(lgrid)
      !common /tbfunc/ tpi,ypi,tpi2,xt0,xt1,xt2,xt3
      !real*8 eav,fsof,plm,qmin,qmax
      !common /pionic/ eav,fsof,plm,qmin,qmax
c
      real*8 vpi(6),wpi(6),wpinn(6),wpind(6),wpidd(6)
      integer*4 :: lp,no,np,npi,npf
      integer*4 :: lt,lst,ltnt,i,k,l
      real*8 :: qscale,upi,fs,fsmpi3,upis,plms,x,y,z,eav1,eav2,dq,qq
      real*8 :: totpi,qqs,ffs,ff,fac,qr,aj0qr,aj2qr,tempnn,tempnd,tempdd
      real*8 :: totnn,totnd,totdd
      character*8 label,psilab(2)
      data psilab/'k-space ','r-space '/
      lt=3-nm
      lst=lt+2
      ltnt=lst+2
      qscale=100.
      if (np.eq.3) then
        upi=138.03/197.33
        fs=.081
      else
        upi=.7
        fs=.0757
      end if
      fsmpi3=fs*upi*197.33/3.
      upis=upi*upi
      plms=plm*plm
      label=psilab(2)
      if (npf.eq.-1) then
        label=psilab(1)
        npf=1
        do 10 i=1,lp
        x=upi*r(i)
        y=plm*r(i)
        z=plm/upi
        ypi(i)=exp(-x)/x-z**2*exp(-y)/x-.5*z*(z**2-1)*(1-2/y)*exp(-y)
        tpi(i)=(1+3/x+3/x**2)*exp(-x)/x-z**2*(1+3/y+3/y**2)*exp(-y)/x
     &   -.5*z*(z**2-1)*(1+1/y)*exp(-y)
   10   continue
      end if
c    *************************
c    compute opep contribution
c    *************************
      x=2.*pi*rho*dr
      vpi(lst)=0.
      vpi(ltnt)=0.
      do 20 i=1,lp
      vpi(lst)=vpi(lst)+ypi(i)*gnn(i,lst)*rs(i)
      vpi(ltnt)=vpi(ltnt)+tpi(i)*gnn(i,ltnt)*rs(i)
   20 continue
      vpi(lst)=fsmpi3*aa(lst)*vpi(lst)*x
      vpi(ltnt)=fsmpi3*aa(ltnt)*vpi(ltnt)*x
      vpi(1)=vpi(lst)+vpi(ltnt)
c print ================================================================
      if (no.ge.1) then
      write(nlog,1495) vpi(1),vpi(lst),vpi(ltnt)
      write(nout,1495) vpi(1),vpi(lst),vpi(ltnt)
 1495 format(/4x,'<ope> =',f8.3,' with <st> =',f8.3,' and <tnt> =',f8.3)
      write(nlog,1500) label,fsof,eav,npf,plm
      write(nout,1500) label,fsof,eav,npf,plm
 1500 format(/4x,'% excess pions/nucleon calculated in closure ',
     & 'approximation',/4x,
     & 'with ',a8,'wave function and (f*/f)**2 =',f6.3,', eav =',f6.1,
     & /4x,'form factor = ( (lam**2-mu**2)/(lam**2+q**2) )**',i1,
     & ' with lam =',f6.3,' fm-1',
     & //4x,'q',5x,'npi(q)',4x,
     & 'c',7x,'t',7x,'s',7x,'st',6x,'tn',6x,'tnt')
      end if
c ======================================================================
      eav1=eav+293.
      eav2=eav1+293.
      dq=(qmax-qmin)/(npi-1)
      qq=qmin-dq
      totpi=0.
      do 100 k=1,npi
      qq=qq+dq
      qqs=qq*qq
      ff=((plms-upis)/(qqs+plms))**npf
      ffs=ff*ff
      fac=4.*fs*qqs*qqs*ffs/(6.*pi*upis*(qqs+upis)**1.5)
      fac=2.*pi*rho*dr*fac*qscale
      do 30 l=1,6,nm
      wpinn(l)=0.
      wpind(l)=0.
   30 wpidd(l)=0.
      do 40 i=1,lp
      qr=qq*r(i)
      aj0qr=sin(qr)/qr
      aj2qr=sin(qr)*(3./qr**3-1./qr)-3.*cos(qr)/qr**2
      x=fsmpi3*fac*(ypi(i)*aj0qr-2.*tpi(i)*aj2qr)*rs(i)
      wpind(1)=wpind(1)-x*(16.*fsof/eav1)*gnn(i,1)
      wpidd(1)=wpidd(1)-x*(32.*fsof**2/(9.*eav2))*gnn(i,1)
      wpind(lt)=wpind(lt)-x*(16.*fsof/(3.*eav1))*gnn(i,lt)*aa(lt)
      wpidd(lt)=wpidd(lt)+x*(16.*fsof**2/(27.*eav2))*gnn(i,lt)*aa(lt)
      x=fsmpi3*fac*(ypi(i)*aj0qr+tpi(i)*aj2qr)*rs(i)
      wpind(3)=wpind(3)-x*(16.*fsof/(3.*eav1))*gnn(i,3)*aa(3)
      wpidd(3)=wpidd(3)+x*(16.*fsof**2/(27.*eav2))*gnn(i,3)*aa(3)
      wpinn(lst)=wpinn(lst)+fac*aj0qr*rs(i)*gnn(i,lst)*aa(lst)
      wpind(lst)=wpind(lst)-x*(16.*fsof/(9.*eav1))*gnn(i,lst)*aa(lst)
      wpidd(lst)=wpidd(lst)-x*(8.*fsof**2/(81.*eav2))*gnn(i,lst)*aa(lst)
      x=fsmpi3*fac*(tpi(i)*aj0qr-ypi(i)*aj2qr+2.*tpi(i)*aj2qr)*rs(i)
      wpind(5)=wpind(5)+x*(8.*fsof/(3.*eav1))*gnn(i,5)*aa(5)
      wpidd(5)=wpidd(5)-x*(8.*fsof**2/(27.*eav2))*gnn(i,5)*aa(5)
      wpinn(ltnt)=wpinn(ltnt)-fac*aj2qr*rs(i)*gnn(i,ltnt)*aa(ltnt)
      wpind(ltnt)=wpind(ltnt)+x*(8.*fsof/(9.*eav1))*gnn(i,ltnt)*aa(ltnt)
      wpidd(ltnt)=wpidd(ltnt)+x*(4.*fsof**2/(81.*eav2))*gnn(i,ltnt)
     & *aa(ltnt)
   40 continue
      temp=0.
      tempnn=0.
      tempnd=0.
      tempdd=0.
      do 50 l=1,6,nm
      wpi(l)=wpinn(l)+wpind(l)+wpidd(l)
      temp=temp+wpi(l)
      tempnn=tempnn+wpinn(l)
      tempnd=tempnd+wpind(l)
   50 tempdd=tempdd+wpidd(l)
      totpi=totpi+temp
      totnn=totnn+tempnn
      totnd=totnd+tempnd
      totdd=totdd+tempdd
c print ================================================================
      if (.not.(no.eq.0)) then
      write(nlog,1505) qq,temp,wpi
      write(nout,1505) qq,temp,wpi
 1505 format(8f8.3)
      if (no.ge.2) then
        write(nout,1510) tempnn,wpinn
        write(nout,1510) tempnd,wpind
        write(nout,1510) tempdd,wpidd
 1510 format(8x,7f8.3)
      end if
      end if
c ======================================================================
  100 continue
      totpi=totpi*dq
      totnn=totnn*dq
      totnd=totnd*dq
      totdd=totdd*dq
c print ================================================================
      if (.not.(no.eq.0)) then
      write(nlog,1515) totpi
      write(nout,1515) totpi
 1515 format(/4x,'total excess = ',f8.3,' % pions/nucleon')
      if (no.eq.2) then
        write(nout,1520) totnn,totnd,totdd
 1520 format(4x,'of which nn part = ',f8.3,' nd part = ',f8.3
     &,' dd part = ',f8.3)
      end if
      end if
c ======================================================================
      return
      end
c *id* nmout ***********************************************************
c subroutine for printing out f, fp, fds, v, g, etc.
c **********************************************************************
      subroutine nmout(le,lg,lt,l3,nie,no,nt,nv)
      use nmvar
      implicit none
      !implicit real*8 (a-h,o-z)
      !implicit integer*4 (i-n)
      !include "params.f"
      !parameter (legrid=lgrid*(lgrid**2+1)/2)
      !parameter (nlog=0,nin=5,nout=6)
      !real*8 kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,h2m,h2mcs,pi,s
      !common /consts/ kf,rho,acn,ast,atn,als,cn,cne,dt,dr,evx,
      !&       h2m,h2mcs,pi,s
      !real*8 r(lgrid),ri(lgrid),rs(lgrid),sl(lgrid),sls(lgrid),
      !&       slp(lgrid),slps(lgrid),sldp(lgrid),sltp(lgrid),
      !&       rllp(lgrid),rlssx(lgrid),rsdsl(lgrid)
      !common /rslate/ r,ri,rs,sl,sls,slp,slps,sldp,sltp,rllp,rlssx,rsdsl
      !real*8 f(lgrid,8),fp(lgrid,8),fds(lgrid,8),v(lgrid,14)
      !common /correl/ f,fp,fds,v
      !real*8 gca(lgrid,6),gcb(lgrid,6),gdd(lgrid,6),gde(lgrid,6),
      !&       gee(lgrid,6),gl(lgrid),gx(lgrid),gy(lgrid),gz(lgrid),
      !&       gnn(lgrid,14)
      !common /gchain/ gca,gcb,gdd,gde,gee,gl,gx,gy,gz,gnn
      !real*8 gfdd(lgrid,6),gfde(lgrid,6),gfed(lgrid,6),gfcc(lgrid,6),
      !&       ghdd(lgrid,6),ghde(lgrid,6),ghed(lgrid,6),ghcc(lgrid,6),
      !&       grdc(lgrid,6),grdd(lgrid,6),grde(lgrid,6),gred(lgrid,6),
      !&       gree(lgrid,6),grfc(lgrid,6),grfd(lgrid,6),grfe(lgrid,6),
      !&       grmd(lgrid,6),grme(lgrid,6)
      !real*8 eca(lgrid,6),ecb(lgrid,6),edd(lgrid,6),ede(lgrid,6),
      !&       eee(lgrid,6),sddd(legrid),sdde(legrid),
      !&       sdee(legrid),seee(legrid),sdcc(legrid),secc(legrid)
      !common /echain/ eca,ecb,edd,ede,eee,sddd,sdde,sdee,seee,sdcc,secc
      !common /mocfun/ gfdd,gfde,gfed,gfcc,ghdd,ghde,ghed,ghcc,
      !&       grdc,grdd,grde,gred,gree,grfc,grfd,grfe,grmd,grme
      !real*8 v3cc(lgrid,6),v3dd(lgrid,6),v3de(lgrid,6),v3ee(lgrid,6)
      !common /tbpots/ v3cc,v3dd,v3de,v3ee
c
      real*8 f00(lgrid),f01(lgrid),f10(lgrid),f11(lgrid),
     &       ft0(lgrid),ft1(lgrid)
      real*8 g01(lgrid),g11(lgrid),h01(lgrid),h11(lgrid)
      integer*4 :: le,lg,lt,l3,nie,no,nt,nv
      integer*4 :: lprt,j,i
      f00(:)=f(:,1)-3*f(:,2)-3*f(:,3)+9*f(:,4)
      f01(:)=f(:,1)+  f(:,2)-3*f(:,3)-3*f(:,4)
      f10(:)=f(:,1)-3*f(:,2)+  f(:,3)-3*f(:,4)
      f11(:)=f(:,1)+  f(:,2)+  f(:,3)+  f(:,4)
      ft0(:)=f(:,5)-3*f(:,6)
      ft1(:)=f(:,5)+  f(:,6)
      g01(:)=gnn(:,1)+gnn(:,2)-3*(gnn(:,3)+gnn(:,4))
      g11(:)=gnn(:,1)+gnn(:,2)+  (gnn(:,3)+gnn(:,4))
      h01(:)=g01(:)/(1+sl(:)**2)
      h11(:)=g11(:)/(1-sl(:)**2)
c
      lprt=lt+1
c   ---------------
c   output f,fp,fds,v
c   ---------------
      write(nout,1998)
 1998 format(/3x,'r',10x,'fST(r)')
      write(nout,1999)
 1999 format(14x,'00',9x,'01',9x,'10',9x,'11',9x,'t0',9x,'t1')
      write(nout,2011) (r(j),f00(j),f01(j),f10(j),f11(j),ft0(j),ft1(j)
     &                 ,j=1,lprt)
      write(nout,2000)
 2000 format(/3x,'r',10x,'f(r,p)')
      write(nout,2006)
 2006 format(14x,'c',10x,'t',10x,'s',10x,'st',9x,'tn',9x,'tnt')
      write(nout,2011) (r(j),(f(j,i),i=1,6),j=1,lprt)
 2011 format(1x,0p,f6.3,5x,1p,6e11.4)
      write(nout,2020)
 2020 format(/3x,'r',10x,'fp(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(fp(j,i),i=1,6),j=1,lprt)
      write(nout,2030)
 2030 format(/3x,'r',10x,'fds(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(fds(j,i),i=1,6),j=1,lprt)
      if (nv.gt.6) then
        write(nout,2040)
 2040   format(/3x,'r',10x,'f(r,p)',16x,'fp(r,p)',15x,'fds(r,p)')
        write(nout,2048)
 2048   format(14x,'b',10x,'bt',9x,'b',10x,'bt',9x,'b',10x,'bt')
        write(nout,2011) (r(j),(f(j,i),i=7,8),(fp(j,i),i=7,8)
     &                   ,(fds(j,i),i=7,8),j=1,lprt)
      end if
      write(nout,2050)
 2050 format(/3x,'r',10x,'v(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(v(j,i),i=1,6),j=1,lprt)
      if (nv.gt.6) then
        write(nout,2050)
        write(nout,2078)
 2078   format(14x,'b',10x,'bt',9x,'q',10x,'qt',9x,'qs',9x,'qst')
        write(nout,2011) (r(j),(v(j,i),i=7,12),j=1,lprt)
      end if
      if (nv.gt.12) then
        write(nout,2050)
        write(nout,2079)
 2079   format(14x,'bb',9x,'bbt')
        write(nout,2082) (r(j),(v(j,i),i=13,14),j=1,lprt)
 2082   format(1x,0p,f6.3,5x,1p,2e11.4)
      end if
      write(nout,2090)
 2090 format(/3x,'r',10x,'gx',9x,'gl',9x,'sl',9x,'slp',8x,'sldp'
     &       ,7x,'sltp')
      write(nout,2011) (r(j),gx(j),gl(j),sl(j),slp(j),sldp(j),sltp(j),
     & j=1,2*lt)
      write(nout,3000)
 3000 format(/3x,'r',10x,'gnn')
      write(nout,2006)
      write(nout,2011) (r(j),gnn(j,1)-1,(gnn(j,i),i=2,6),j=1,2*lt)
      if (nv.gt.6) then
        write(nout,3000)
        write(nout,2078)
        write(nout,2011) (r(j),(gnn(j,i),i=7,12),j=1,2*lt)
      end if
      if (nv.gt.12) then
        write(nout,3000)
        write(nout,2079)
        write(nout,2082) (r(j),(gnn(j,i),i=13,14),j=1,2*lt)
      end if
c   ---------------
      if (.not.(no.le.2)) then
      write(nout,3001)
 3001 format(/3x,'r',10x,'f01**2',5x,'f11**2',5x,'g01',8x,'g11',8x
     &                  ,'h01',8x,'h11')
      write(nout,2011) (r(j),f01(j)**2,f11(j)**2,g01(j),g11(j)
     &                      ,h01(j),h11(j),j=1,2*lt)
c   ---------------
c   output g chains
c   ---------------
      if (.not.(no.le.3)) then
      write(nout,3005)
 3005 format(/3x,'r',10x,'gdd(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(gdd(j,i),i=1,6),j=1,lg)
      write(nout,3010)
 3010 format(/3x,'r',10x,'gde(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(gde(j,i),i=1,6),j=1,lg)
      write(nout,3020)
 3020 format(/3x,'r',10x,'gee(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(gee(j,i),i=1,6),j=1,lg)
      write(nout,3030)
 3030 format(/3x,'r',10x,'gca(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(gca(j,i),i=1,6),j=1,lg)
      write(nout,3040)
 3040 format(/3x,'r',10x,'gcb(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(gcb(j,i),i=1,6),j=1,lg)
      if (nie.ge.1) then
        write(nout,3050)
 3050   format(/3x,'r',10x,'edd(r,p)')
        write(nout,2006)
        write(nout,2011) (r(j),(edd(j,i),i=1,6),j=1,le)
        write(nout,3060)
 3060   format(/3x,'r',10x,'ede(r,p)')
        write(nout,2006)
        write(nout,2011) (r(j),(ede(j,i),i=1,6),j=1,le)
        write(nout,3070)
 3070   format(/3x,'r',10x,'eee(r,p)')
        write(nout,2006)
        write(nout,2011) (r(j),(eee(j,i),i=1,6),j=1,le)
        write(nout,3080)
 3080   format(/3x,'r',10x,'eca(r,p)')
        write(nout,2006)
        write(nout,2011) (r(j),(eca(j,i),i=1,6),j=1,le)
        write(nout,3090)
 3090   format(/3x,'r',10x,'ecb(r,p)')
        write(nout,2006)
        write(nout,2011) (r(j),(ecb(j,i),i=1,6),j=1,le)
      end if
      write(nout,3105)
 3105 format(/3x,'r',10x,'grdd(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(grdd(j,i),i=1,6),j=1,lg)
      write(nout,3110)
 3110 format(/3x,'r',10x,'grde(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(grde(j,i),i=1,6),j=1,lg)
      write(nout,3115)
 3115 format(/3x,'r',10x,'grfd(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(grfd(j,i),i=1,6),j=1,lg)
      write(nout,3120)
 3120 format(/3x,'r',10x,'grfe(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(grfe(j,i),i=1,6),j=1,lg)
      write(nout,3205)
 3205 format(/3x,'r',10x,'gfdd(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(gfdd(j,i),i=1,6),j=1,lg)
      write(nout,3210)
 3210 format(/3x,'r',10x,'gfde(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(gfde(j,i),i=1,6),j=1,lg)
      write(nout,3215)
 3215 format(/3x,'r',10x,'gfed(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(gfed(j,i),i=1,6),j=1,lg)
      write(nout,3220)
 3220 format(/3x,'r',10x,'gfcc(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(gfcc(j,i),i=1,6),j=1,lg)
      write(nout,3305)
 3305 format(/3x,'r',10x,'ghdd(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(ghdd(j,i),i=1,6),j=1,lg)
      write(nout,3310)
 3310 format(/3x,'r',10x,'ghde(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(ghde(j,i),i=1,6),j=1,lg)
      write(nout,3315)
 3315 format(/3x,'r',10x,'ghed(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(ghed(j,i),i=1,6),j=1,lg)
      write(nout,3320)
 3320 format(/3x,'r',10x,'ghcc(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(ghcc(j,i),i=1,6),j=1,lg)
c   -----------------
c   ouput v3xx chains
c   -----------------
   10 if (.not.(nt.eq.0.or.nt.ge.4)) then
      write(nout,5000)
 5000 format(/3x,'r',10x,'v3dd(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(v3dd(j,i,1),i=1,6),j=1,l3)
      write(nout,5010)
 5010 format(/3x,'r',10x,'v3de(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(v3de(j,i,1),i=1,6),j=1,l3)
      write(nout,5020)
 5020 format(/3x,'r',10x,'v3ee(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(v3ee(j,i,1),i=1,6),j=1,l3)
      write(nout,5030)
 5030 format(/3x,'r',10x,'v3cc(r,p)')
      write(nout,2006)
      write(nout,2011) (r(j),(v3cc(j,i,1),i=1,6),j=1,l3)
      end if
      end if
      end if
      return
      end
