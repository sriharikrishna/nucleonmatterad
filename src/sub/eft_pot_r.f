      subroutine setpot_chiral(lpot,h2m,h2mcsb)
c
      implicit real*8(a-h,o-z)
      save
c
      include "params_new_pot.f"
      parameter( nrl=nr/2 )
c
      common/gf/gx(nz),wgx(nz)
      common/gp/gg(nx),wgg(nx)
c
      common/lecs/rcf,c0c0,c0c1,c2c00,c2c10,c2c01,c2c11,
     x                           c4c00,c4c10,c4c01,c4c11,
     x                c2t0,c2t1,c4t0,c4t1,c2b,c4b0,c4b1,c4bb,
     x                c4q0,c4q1,c4p0,c4p1,c4tp0,c4tp1,
     x                c0cv,c2c0v,c2c1v,c2tv,c2bv,
     x                c0ct,c2c0t,c2c1t,c2tt,c2bt,
     x                ccc1,ccc2,ccc3,ccc4,cb38,chan0
      common/lecs_st/cc0cx  (0:1,0:1),cc2cx  (0:1,0:1),cc4cx (0:1,0:1),
     x               ct2tx  (0:1,0:1),ct4tx  (0:1,0:1),cb2bx (0:1,0:1),
     x               cb4bx  (0:1,0:1),cq4qx  (0:1,0:1),cp4px (0:1,0:1),
     x               cbb4bbx(0:1,0:1),ctp4tpx(0:1,0:1)
c
      common/lecs_st_cd/cc0cvx(0:1,0:1),cc2cvx(0:1,0:1),ct2tvx(0:1,0:1),
     x                                                  cb2bvx(0:1,0:1),
     x                  cc0ctx(0:1,0:1),cc2ctx(0:1,0:1),ct2ttx(0:1,0:1),
     x                                                  cb2btx(0:1,0:1)
c
      logical lpotls,lpotll
      common /logpot/ lpotls,lpotll
c
      dimension vnlo(8), vdnlo(6,2), vn2lo(6),vdn2lo(6,2)
      dimension cc0cvxpp(0:1,0:1),cc2cvxpp(0:1,0:1),ct2tvxpp(0:1,0:1),
     x          cb2bvxpp(0:1,0:1),cc0ctxpp(0:1,0:1),cc2ctxpp(0:1,0:1),
     x                            ct2ttxpp(0:1,0:1),cb2btxpp(0:1,0:1)
      dimension cc0ctxnp(0:1,0:1),cc2ctxnp(0:1,0:1),ct2ttxnp(0:1,0:1),
     x                                              cb2btxnp(0:1,0:1)
      dimension cc0cvxnn(0:1,0:1),cc2cvxnn(0:1,0:1),ct2tvxnn(0:1,0:1),
     x          cb2bvxnn(0:1,0:1),cc0ctxnn(0:1,0:1),cc2ctxnn(0:1,0:1),
     x                            ct2ttxnn(0:1,0:1),cb2btxnn(0:1,0:1)
      dimension vv(22),vp(12)
c
      data idlt / 1 /
      data amnp / 4.7549097d0 / , amnn / 4.7614639d0 /

c
      h2m=0.5d0*hbarc*( 1.d0/amnp+1.d0/amnn )
      h2mcsb=0.5d0*hbarc*( 1.d0/amnp-1.d0/amnn )
      lpotls=.true.
      lpotll=.true.
c
      if (lpot .eq. 101) ilb=1
      if (lpot .eq. 102) ilb=2
      if (lpot .eq. 103) ilb=3
      if (lpot .eq. 104) ilb=4
      if (lpot .eq. 105) ilb=5
      if (lpot .eq. 106) ilb=6
      if (lpot .eq. 107) ilb=7
      if (lpot .eq. 108) ilb=8
      if (lpot .eq. 109) ilb=9
      if (lpot .eq. 110) ilb=10
      if (lpot .eq. 111) ilb=11
      if (lpot .eq. 112) ilb=12
      if (lpot .eq. 113) ilb=13
c
      call gauss_pot( nz , 0.d0 , 1.d0 , gx , wgx )
      call gauss_pot( nx , 0.d0 , 1.d0 , gg , wgg )
c
      call set_eft_pot( 1, ilb )
c
      do 10 it=0,1
      do 10 is=0,1
      cc0cvxpp(is,it)=cc0cvx(is,it)
      cc2cvxpp(is,it)=cc2cvx(is,it)
      ct2tvxpp(is,it)=ct2tvx(is,it)
      cb2bvxpp(is,it)=cb2bvx(is,it)
c
      cc0ctxpp(is,it)=cc0ctx(is,it)
      cc2ctxpp(is,it)=cc2ctx(is,it)
      ct2ttxpp(is,it)=ct2ttx(is,it)
      cb2btxpp(is,it)=cb2btx(is,it)
 10   continue
c
      call set_eft_pot( 0, ilb )
c
      do 20 it=0,1
      do 20 is=0,1
      cc0ctxnp(is,it)=cc0ctx(is,it)
      cc2ctxnp(is,it)=cc2ctx(is,it)
      ct2ttxnp(is,it)=ct2ttx(is,it)
      cb2btxnp(is,it)=cb2btx(is,it)
 20   continue
c
      call set_eft_pot(-1, ilb )
c
      do 30 it=0,1
      do 30 is=0,1
      cc0cvxnn(is,it)=cc0cvx(is,it)
      cc2cvxnn(is,it)=cc2cvx(is,it)
      ct2tvxnn(is,it)=ct2tvx(is,it)
      cb2bvxnn(is,it)=cb2bvx(is,it)
c
      cc0ctxnn(is,it)=cc0ctx(is,it)
      cc2ctxnn(is,it)=cc2ctx(is,it)
      ct2ttxnn(is,it)=ct2ttx(is,it)
      cb2btxnn(is,it)=cb2btx(is,it)
 30   continue
c
      return
c************************************************************
      entry pot_chiral(rr,vv,vp)
c
      do i=1,22
      vv(i)=0.d0
      enddo
      
      call chiral_nlo ( idlt , rr , vnlo , vdnlo  )
      call chiral_n2lo( idlt , rr , vn2lo ,vdn2lo )
c
      if ( ilb .eq. 13) then
          vnlo(1)=0.d0
          vnlo(2)=0.d0
          vnlo(3)=0.d0
          vnlo(5)=0.d0
c
         do j=1,6
          vdnlo(j,1)=0.d0
          vdnlo(j,2)=0.d0
         enddo
c
      endif
c
      if( ilb .eq. 12 .or. ilb .eq. 13 ) then
c
         do j=1,6
          vn2lo (j  )=0.d0
          vdn2lo(j,1)=0.d0
          vdn2lo(j,2)=0.d0
         enddo
c
      endif
c
      vclr =vnlo(1)+vdnlo(1,1)+vdnlo(1,2)+vn2lo(1)+vdn2lo(1,1)
     x                                            +vdn2lo(1,2)
      vctlr=vnlo(2)+vdnlo(2,1)+vdnlo(2,2)+vn2lo(2)+vdn2lo(2,1)
     x                                            +vdn2lo(2,2)
      vcslr=vnlo(3)+vdnlo(3,1)+vdnlo(3,2)+vn2lo(3)+vdn2lo(3,1)
     x                                            +vdn2lo(3,2)
      vtslr=vnlo(4)+vdnlo(4,1)+vdnlo(4,2)+vn2lo(4)+vdn2lo(4,1)
     x                                            +vdn2lo(4,2)
c
      vtlr =vnlo(5)+vdnlo(5,1)+vdnlo(5,2)+vn2lo(5)+vdn2lo(5,1)
     x                                            +vdn2lo(5,2)
      vttlr=vnlo(6)+vdnlo(6,1)+vdnlo(6,2)+vn2lo(6)+vdn2lo(6,1)
     x                                            +vdn2lo(6,2)
c
c central part of the long-range part of the potential projected on S,T
c
      vclr00=vclr-3.d0*vctlr-3.d0*vcslr+9.d0*vtslr
      vclr01=vclr+     vctlr-3.d0*vcslr-3.d0*vtslr
      vclr10=vclr-3.d0*vctlr+     vcslr-3.d0*vtslr
      vclr11=vclr+     vctlr+     vcslr+     vtslr
c
c tensor part of the long-range part of the potential projected on S=1,T
c
      vtlr0=vtlr-3.d0*vttlr
      vtlr1=vtlr+     vttlr
c
c isospin symmetry breaking terms of the long-range part of the
c potential projected on S,T=1,Tz
c
      vctlr0pp=-6.d0*vnlo(7)
      vctlr0np=12.d0*vnlo(7)
      vctlr0nn=-6.d0*vnlo(7)
c
      vctlr1pp= 2.d0*vnlo(7)
      vctlr1np=-4.d0*vnlo(7)
      vctlr1nn= 2.d0*vnlo(7)
c
      vttlrpp= 2.d0*vnlo(8)
      vttlrnp=-4.d0*vnlo(8)
      vttlrnn= 2.d0*vnlo(8)
c
c this call defines the radial functions for the short range part of the
c potential
c
      call fnt_cntc( ncf , ilb , rr , fc0  ,
     x                                fc2  , ft2   , fb2 ,
     x                                fc4  , ft4   , fb4 ,
     x                                fbb4 , fq4   ,
     x                                fp4  , fp4d  ,
     x                                ftp4 , ftp4d )
c 
c Here we define the isospin-independent part of the short-range
c potential projected on S,T
c
      vcsr00=cc0cx(0,0)*fc0+cc2cx(0,0)*fc2+cc4cx(0,0)*fc4
      vcsr01=cc0cx(0,1)*fc0+cc2cx(0,1)*fc2+cc4cx(0,1)*fc4
      vcsr10=cc0cx(1,0)*fc0+cc2cx(1,0)*fc2+cc4cx(1,0)*fc4
      vcsr11=cc0cx(1,1)*fc0+cc2cx(1,1)*fc2+cc4cx(1,1)*fc4
c
      vtsr0=ct2tx(1,0)*ft2+ct4tx(1,0)*ft4
      vtsr1=ct2tx(1,1)*ft2+ct4tx(1,1)*ft4
c
      vbsr0=cb2bx(1,0)*fb2+cb4bx(1,0)*fb4
      vbsr1=cb2bx(1,1)*fb2+cb4bx(1,1)*fb4
c
      vqsr0=cq4qx(0,1)*fq4
      vqsr1=cq4qx(1,1)*fq4
c
      vbbsr=cbb4bbx(1,1)*fbb4
c
c ATTENTION: define here the p2 potential
c
      if( ilb .le. 3 ) then
         vpsr0=cp4px(0,1)*fp4
         vpsr1=cp4px(1,1)*fp4
c
         vpdsr0=cp4px(0,1)*fp4d
         vpdsr1=cp4px(1,1)*fp4d
c
         vtpsr0=ctp4tpx(1,0)*ftp4
         vtpsr1=ctp4tpx(1,1)*ftp4
c
         vtpdsr0=ctp4tpx(1,0)*ftp4d
         vtpdsr1=ctp4tpx(1,1)*ftp4d
      endif
c
c Here we define the isospin-dependent part of the short-range potential
c on S,T for pp
c
      vcvsr0pp=cc0cvxpp(0,1)*fc0+cc2cvxpp(0,1)*fc2
      vcvsr1pp=cc0cvxpp(1,1)*fc0+cc2cvxpp(1,1)*fc2
c
      vtvsrpp=ct2tvxpp(1,1)*ft2
c
      vbvsrpp=cb2bvxpp(1,1)*fb2
c
      vctsr0pp=cc0ctxpp(0,1)*fc0+cc2ctxpp(0,1)*fc2
      vctsr1pp=cc0ctxpp(1,1)*fc0+cc2ctxpp(1,1)*fc2
c
      vttsrpp=ct2ttxpp(1,1)*ft2
c
      vbtsrpp=cb2btxpp(1,1)*fb2
c
c Here we define the isospin-dependent part of the short-range potential
c on S,T for np
c
      vctsr0np=cc0ctxnp(0,1)*fc0+cc2ctxnp(0,1)*fc2
      vctsr1np=cc0ctxnp(1,1)*fc0+cc2ctxnp(1,1)*fc2
c      
      vttsrnp=ct2ttxnp(1,1)*ft2
c
      vbtsrnp=cb2btxnp(1,1)*fb2
c
c Here we define the isospin-dependent part of the short-range potential
c on S,T for nn
c
      vcvsr0nn=cc0cvxnn(0,1)*fc0+cc2cvxnn(0,1)*fc2
      vcvsr1nn=cc0cvxnn(1,1)*fc0+cc2cvxnn(1,1)*fc2
c      
      vtvsrnn=ct2tvxnn(1,1)*ft2
c      
      vbvsrnn=cb2bvxnn(1,1)*fb2
c      
      vctsr0nn=cc0ctxnn(0,1)*fc0+cc2ctxnn(0,1)*fc2
      vctsr1nn=cc0ctxnn(1,1)*fc0+cc2ctxnn(1,1)*fc2
c      
      vttsrnn=ct2ttxnn(1,1)*ft2
c
      vbtsrnn=cb2btxnn(1,1)*fb2
c
c Define here the components of the full potential
c
      p00  =vclr00+vcsr00
      p10  =vclr10+vcsr10
      p01pp=vclr01+vctlr0pp+vcsr01+vcvsr0pp+vctsr0pp
      p01np=vclr01+vctlr0np+vcsr01+         vctsr0np
      p01nn=vclr01+vctlr0nn+vcsr01+vcvsr0nn+vctsr0nn
      p11pp=vclr11+vctlr1pp+vcsr11+vcvsr1pp+vctsr1pp
      p11np=vclr11+vctlr1np+vcsr11+         vctsr1np
      p11nn=vclr11+vctlr1nn+vcsr11+vcvsr1nn+vctsr1nn  
c
      pt0  =vtlr0+vtsr0
      pt1pp=vtlr1+vttlrpp+vtsr1+vtvsrpp+vttsrpp
      pt1np=vtlr1+vttlrnp+vtsr1+        vttsrnp
      pt1nn=vtlr1+vttlrnn+vtsr1+vtvsrnn+vttsrnn
c
      pb0  =vbsr0
      pb1pp=vbsr1+vbvsrpp+vbtsrpp
      pb1np=vbsr1+        vbtsrnp
      pb1nn=vbsr1+vbvsrnn+vbtsrnn
c
      pq0=vqsr0
      pq1=vqsr1
c
      pbb=vbbsr
c
c ATTENTION: p2 potential here
c
      if(ilb .le. 3) then
c
       pp0=vpsr0
       pp1=vpsr1
c
       ptp0=vtpsr0
       ptp1=vtpsr1
c
      endif
c       
      p01  =(p01pp+p01np+p01nn)/3.d0
      p01cd=(.5d0*(p01pp+p01nn)-p01np)/6.d0
      p01cs=(p01pp-p01nn)/4.d0
c
      p11  =(p11pp+p11np+p11nn)/3.d0
      p11cd=(.5d0*(p11pp+p11nn)-p11np)/6.d0
      p11cs=(p11pp-p11nn)/4.d0
c
      pt1  =(pt1pp+pt1np+pt1nn)/3.d0
      pt1cd=(.5d0*(pt1pp+pt1nn)-pt1np)/6.d0
      pt1cs=(pt1pp-pt1nn)/4.d0
c
      pb1  =(pb1pp+pb1np+pb1nn)/3.d0
      pb1cd=(.5d0*(pb1pp+pb1nn)-pb1np)/6.d0
      pb1cs=(pb1pp-pb1nn)/4.d0
c
      vv(1)=.0625d0*(9.d0*p11+3.d0*p10+3.d0*p01+p00)
      vv(2)=.0625d0*(3.d0*p11-3.d0*p10+     p01-p00)
      vv(3)=.0625d0*(3.d0*p11+     p10-3.d0*p01-p00)
      vv(4)=.0625d0*(     p11-     p10-     p01+p00)  
c
      vv(5)=.25d0*(3*pt1+pt0)
      vv(6)=.25d0*(  pt1-pt0)
c
      vv(7)=.25d0*(3*pb1+pb0)
      vv(8)=.25d0*(  pb1-pb0)
c
      vv(9) =.25d0*(3*pq1+pq0)
      vv(10)=.0d0
      vv(11)=.25d0*(  pq1-pq0)
      vv(12)=.0d0
c
      vv(13)=pbb
      vv(14)=0.d0
c
      vv(15)=.25d0*(3*p11cd+p01cd)
      vv(16)=.25d0*(  p11cd-p01cd)
c
      vv(17)=pt1cd
c
      vv(18)=pb1cd        
c
      if (ilb .ge. 4) then
      vv(19)=p01cs
      else
      vv(19)=.25d0*(3*p11cs+p01cs)
      vv(20)=.25d0*(  p11cs-p01cs)
c
      vv(21)=pt1cs
c
      vv(22)=pb1cs
      endif


 40   continue
c
      return
      end
c_______________________________________________________________________
      subroutine set_eft_pot( itz, ilb)
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      character*2 lm
      character*5 clecs
      save
c
      include "params_new_pot.f"
      parameter ( mpx=41 )
c
      common/lecs/rcf,c0c0,c0c1,c2c00,c2c10,c2c01,c2c11,
     x                          c4c00,c4c10,c4c01,c4c11,
     x                c2t0,c2t1,c4t0,c4t1,c2b,c4b0,c4b1,c4bb,
     x                c4q0,c4q1,c4p0,c4p1,c4tp0,c4tp1,
     x                c0cv,c2c0v,c2c1v,c2tv,c2bv,
     x                c0ct,c2c0t,c2c1t,c2tt,c2bt,
     x                ccc1,ccc2,ccc3,ccc4,cb38,chan0
      common/lecs_st/cc0cx  (0:1,0:1),cc2cx  (0:1,0:1),cc4cx (0:1,0:1),
     x               ct2tx  (0:1,0:1),ct4tx  (0:1,0:1),cb2bx (0:1,0:1),
     x               cb4bx  (0:1,0:1),cq4qx  (0:1,0:1),cp4px (0:1,0:1),
     x               cbb4bbx(0:1,0:1),ctp4tpx(0:1,0:1)
c
      common/lecs_st_cd/cc0cvx(0:1,0:1),cc2cvx(0:1,0:1),ct2tvx(0:1,0:1),
     x                                                  cb2bvx(0:1,0:1),
     x                  cc0ctx(0:1,0:1),cc2ctx(0:1,0:1),ct2ttx(0:1,0:1),
     x                                                  cb2btx(0:1,0:1)
c
      dimension lm(mpx),dd(mpx),clecs(mpx),ddr(mpx,13)
      equivalence ( dd( 1),rcf ),
     x            ( dd( 2),c0c0  ),( dd( 3),c0c1  ),
     x            ( dd( 4),c2c00 ),( dd( 5),c2c10 ),
     x            ( dd( 6),c2c01 ),( dd( 7),c2c11 ),
     x            ( dd( 8),c4c00 ),( dd( 9),c4c10 ),
     x            ( dd(10),c4c01 ),( dd(11),c4c11 ),
     x            ( dd(12),c2t0  ),( dd(13),c2t1  ),
     x            ( dd(14),c4t0  ),( dd(15),c4t1  ),
     x            ( dd(16),c2b   ),( dd(17),c4b0  ),
     x            ( dd(18),c4b1  ),( dd(19),c4bb  ),
     x            ( dd(20),c4q0  ),( dd(21),c4q1  ),
     x            ( dd(22),c4p0  ),( dd(23),c4p1  ),
     x            ( dd(24),c4tp0 ),( dd(25),c4tp1 ),
     x            ( dd(26),c0cv  ),( dd(27),c2c0v ),
     x            ( dd(28),c2c1v ),( dd(29),c2tv  ),
     x            ( dd(30),c2bv  ),( dd(31),c0ct  ),
     x            ( dd(32),c2c0t ),( dd(33),c2c1t ),
     x            ( dd(34),c2tt  ),( dd(35),c2bt  ),
     x            ( dd(36),ccc1  ),( dd(37),ccc2  ),
     x            ( dd(38),ccc3  ),( dd(39),ccc4  ),
     x            ( dd(40),cb38  ),( dd(41),chan0 )
c
c      read (5,5000) ( lm(i),dd(i),clecs(i),i=1,mpx )

c
      data ddr/0.800000000E+00_8,
     x         0.533703466E+02_8, 0.167269044E+02_8, 0.223074210E+01_8,
     x         0.488868534E+00_8, 0.544889452E+00_8, -.593238340E+00_8,
     x         0.345909595E-01_8, -.653724085E-01_8, 0.231236234E-01_8,
     x                                               -.353649907E-01_8,
     x         0.777038048E-01_8, -.327767098E-01_8, 0.148589889E-01_8,
     x                                               -.215104370E-02_8,
     x         0.756873223E+00_8, 0.761517957E-01_8, 0.312481509E-01_8,
     x                                               0.496626324E-01_8,
     x         -.933318621E-01_8, -.185978150E+00_8, 
     x         0.483841490E-01_8, 0.553625548E-01_8,                                              
     x         -.171856001E-02_8, 0.101982280E-01_8,
     x         0.619532364E+00_8, 0.460824232E-01_8, 0.136321035E-01_8,
     x         -.611590202E-01_8, -.153321183E+01_8,
     x         0.702062973E-01_8, 0.400684414E-02_8, 0.882903390E-02_8,
     x         0.210971577E-01_8, 0.495595154E+00_8,
     x         -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8, 
     x         0.133000000E+01_8, 0.140000000E+01_8,
     x         0.273650324E+01_8, 
c
     x        0.100000000E+01_8, 
     x        0.213459770E+02_8, 0.467382588E+01_8, 0.121844815E+01_8, 
     x        0.477247353E+00_8, 0.504551482E+00_8, -.629190318E+00_8,
     x        -.616841413E-01_8, -.732674011E-01_8, -.471711650E-01_8,
     x                                              -.772080850E-01_8,
     x        0.372927574E-01_8, -.748980434E-01_8, 0.105314022E-01_8,
     x                                              -.499759063E-03_8,
     x        0.733861106E+00_8, 0.101454626E+00_8, -.144601338E-01_8,
     x                                              0.390907294E-01_8,
     x        -.217113694E+00_8, -.202443182E+00_8, 
     x        0.132065078E+00_8, 0.925031384E-01_8,
     x        -.421404352E-02_8, 0.549304077E-02_8,
     x        -.116219208E+00_8, -.177661482E-01_8, 0.737686633E-01_8,
     x        -.770713070E-02_8, -.158113670E+01_8,
     x        0.666916664E-01_8, 0.668128491E-02_8, -.338677404E-01_8,
     x        -.109811363E-03_8, 0.518036830E+00_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8, 
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,   
c
     x        0.120000000E+01_8,
     x        0.698590000E+01_8, 0.342929719E+00_8, 0.472724849E+00_8,
     x        0.419425560E+00_8, 0.350574156E+00_8, -.622004743E+00_8,
     x        -.829833997E-01_8, -.703849975E-01_8, -.729738770E-01_8,
     x                                              -.111853393E+00_8,
     x        0.875966683E-01_8, -.146717246E+00_8, 0.226561861E-01_8,
     x                                              -.396266139E-02_8,
     x        -.114092338E+00_8, 0.638143052E-01_8, -.521790258E-01_8,
     x                                              0.263439181E-02_8,
     x        -.231466008E+00_8, -.161114625E+00_8,
     x        0.162150108E+00_8, 0.961540375E-01_8,
     x        -.116324415E-01_8, 0.653063138E-02_8,
     x        -.873029912E-01_8, -.356288662E-01_8, 0.104690576E+00_8,
     x        0.386207692E-01_8, -.761783605E+00_8,
     x        0.580466231E-01_8, 0.159406714E-01_8, -.370798431E-01_8,
     x        -.139937121E-01_8, 0.258260695E+00_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8, 
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,
c
     x        0.100000000E+01_8,
     x        0.116493463E+02_8, 0.683438283E+01_8, 0.679110717E+00_8,
     x        0.369528990E+00_8, -.370119939E+00_8, -.486045732E+00_8,
     x        0.638354629E-01_8, -.508524690E-01_8, -.816884149E-02_8,
     x                                              -.514106129E-01_8,
     x        0.203490986E+00_8, -.881264853E-01_8, 0.269304858E-01_8,
     x                                              0.332104549E-02_8,
     x        -.119745217E+01_8, -.596572113E-01_8, 0.804170579E-02_8,
     x                                              0.309485137E-01_8,
     x        -.458771185E-01_8, -.923595131E-01_8, 
     x        0.000000000E+00_8, 0.000000000E+00_8,                            
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.101898890E-01_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.241659108E-01_8, 0.709830242E-02_8, -.101019419E-01_8,
     x        -.357044584E-02_8, 0.116703396E-01_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8, 
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,
c
     x        0.100000000E+01_8,
     x        0.116493463E+02_8, 0.683438282E+01_8, 0.679110717E+00_8,
     x        0.369528990E+00_8, -.370119939E+00_8, -.480668142E+00_8,
     x        0.638354629E-01_8, -.508524690E-01_8, -.816884149E-02_8,
     x                                              -.553472136E-01_8,
     x        0.203490986E+00_8, -.927733547E-01_8, 0.269304858E-01_8, 
     x                                              0.298550753E-02_8,
     x        -.119745217E+01_8, -.596572113E-01_8, 0.116454740E-01_8,
     x                                              0.309485137E-01_8,
     x        -.458771185E-01_8, -.923595131E-01_8, 
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.101898890E-01_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.241659108E-01_8, 0.709830242E-02_8, -.730929539E-02_8,
     x        -.455900560E-02_8, 0.185999731E-01_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8,
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,
c
     x        0.120000000E+01_8,
     x        0.429467518E+01_8, 0.220329600E+01_8, 0.208392044E+00_8,
     x        0.435703473E+00_8, -.626058088E-01_8, -.271837681E+00_8,
     x        0.336791587E-01_8, -.456372905E-01_8, -.124476019E-02_8,
     x                                              -.636916547E-01_8,
     x        0.220117703E+00_8, -.820792453E-01_8, 0.294199793E-01_8,
     x                                              -.262822225E-02_8,
     x        -.101720571E+01_8, -.645879980E-01_8, -.926256154E-02_8,
     x                                              0.479701222E-01_8,
     x        -.720887094E-01_8, -.130220515E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.932547676E-02_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.157823996E-01_8, 0.166861667E-01_8, -.846799203E-02_8,
     x        -.579980308E-02_8, 0.225016667E-01_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8,
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,
c
     x        0.800000000E+00_8,
     x        0.369387786E+02_8, 0.124651543E+02_8, 0.146843747E+01_8,
     x        0.172505853E+00_8, -.349240562E+00_8, -.110411408E+01_8,
     x        0.234127471E-01_8, -.481122978E-01_8, 0.222444631E-01_8,
     x                                              -.324358424E-01_8,
     x        0.167255147E+00_8, -.941752141E-01_8, 0.142182983E-01_8,
     x                                              0.523011253E-02_8,
     x        -.146474755E+01_8, -.509326050E-01_8, 0.652002036E-01_8,
     x                                              0.917591033E-01_8,
     x        -.428215305E-01_8, -.149684562E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.135781790E-01_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.219588054E-01_8, 0.116624694E-02_8, -.398644723E-02_8,
     x        -.312608883E-03_8, 0.898753834E-02_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8,
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,
c
     x        0.800000000E+00_8,
     x        0.369325803E+02_8, 0.124668847E+02_8, 0.154851838E+01_8,
     x        0.171531460E+00_8, -.351790507E+00_8, -.110139061E+01_8,
     x        0.416690410E-01_8, -.480226344E-01_8, 0.223603910E-01_8,
     x                                              -.341917808E-01_8,
     x        0.169358784E+00_8, -.919385126E-01_8, 0.142087627E-01_8,
     x                                              0.491451535E-02_8,
     x        -.146587493E+01_8, -.580711243E-01_8, 0.602283869E-01_8,
     x                                              0.913519424E-01_8,
     x        -.276793955E-01_8, -.150302175E+00_8, 
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.131271163E-01_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.139472303E-01_8, 0.271855483E-03_8, -.120442146E-01_8,
     x        0.390586653E-04_8, 0.884404332E-03_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8,
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,
c
     x        0.100000000E+01_8,
     x        0.120215871E+02_8, 0.719080298E+01_8, 0.582589536E+00_8,
     x        0.339139093E+00_8, -.365653112E+00_8, -.469123203E+00_8,
     x        0.658293165E-01_8, -.200563447E-01_8, 0.313476208E-03_8,
     x                                              -.219529091E-01_8,
     x        0.130134408E+00_8, -.915688512E-01_8, 0.298978917E-01_8,
     x                                              0.446309944E-03_8,
     x        -.101984892E+01_8, -.368880746E-01_8, 0.675136686E-02_8,
     x                                              0.609585814E-03_8,
     x        -.387621777E-01_8, -.595078638E-01_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.107754067E-01_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.210214026E-01_8, 0.429062717E-02_8, -.127651664E-02_8,
     x        -.319445941E-03_8, 0.287987270E-02_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8,
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,
c
     x        0.120000000E+01_8,
     x        0.441620986E+01_8, 0.244265097E+01_8, 0.294761748E+00_8,
     x        0.344555359E+00_8, -.874768018E-01_8, -.298615166E+00_8,
     x        0.164469617E-01_8, -.290985295E-01_8, 0.719353106E-02_8,
     x                                              -.496870949E-01_8,
     x        0.178646941E+00_8, -.684264475E-01_8, 0.368976476E-01_8,
     x                                              0.250222042E-02_8,
     x        -.957819135E+00_8, -.147256317E+00_8, -.224866705E-01_8, 
     x                                              0.328163634E-01_8,
     x        -.514376361E-01_8, -.981491616E-01_8, 
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.957569475E-02_8, 0.000000000E+00_8, 0.000000000E+00,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.219475836E-01_8, 0.235135370E-01_8, -.990517998E-02,
     x        -.668274557E-02_8, 0.127697053E-01_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8,
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,
c
     x        0.100000000E+01_8,
     x        0.120047020E+02_8, 0.764367866E+01_8, 0.973877619E+00_8,
     x        -.405996279E+00_8, -.367245385E+00_8, -.757145843E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x                                              0.000000000E+00_8,
     x        0.371544187E+00_8, -.621931538E-01_8, 0.000000000E+00_8,
     x                                              0.000000000E+00_8,
     x        -.750197104E+00_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x                                              0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.102009602E-01_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.277531796E-01_8, 0.694004508E-02_8, -.366844091E-01_8,
     x        -.276363587E-01_8, 0.880784030E-02_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8,
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,
c
     x        0.100000000E+01_8,
     x        0.854300448E+01_8, 0.716183663E+01_8, 0.520823756E+00_8,
     x        -.154705018E+00_8, -.121180643E+00_8, -.534240332E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8, 0.000000000E+00_8, 
     x                                              0.000000000E+00_8,
     x        0.367928965E+00_8, -.266896895E-01_8, 0.000000000E+00_8,
     x                                              0.000000000E+00_8,
     x        -.903141771E+00_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x                                              0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8, 
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.107010233E-01_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.159975754E-02_8, -.331512687E-02_8, -.348957317E-01_8,
     x        -.299812806E-01_8, 0.173082404E-01_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8,
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8,
c
     x        0.100000000E+01_8,
     x        -.184815192E+01_8, -.122857613E+01_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x                                              0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x                                              0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x                                              0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        0.135000070E-02_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        -.303864937E-02_8, 0.000000000E+00_8, 0.000000000E+00_8,
     x        0.000000000E+00_8, 0.000000000E+00_8,
     x        -.570000000E+00_8, -.250000000E+00_8, -.790000000E+00_8,
     x        0.133000000E+01_8, 0.140000000E+01_8,
     x        0.273650324E+01_8/


      do i=1,mpx
      dd(i)=ddr(i,ilb)
      enddo
c
      write(6,5010) ( dd(i),i=1,mpx )

c
c   S=0 and T=0 channel LECs follow
c
      cc0cx(0,0)=c0c0
      cc2cx(0,0)=c2c00
      cc4cx(0,0)=c4c00
c
      ct2tx(0,0)=0.d0
      ct4tx(0,0)=0.d0
c
      cb2bx(0,0)=0.d0
      cb4bx(0,0)=0.d0
c
      cq4qx(0,0)=c4q0
      cp4px(0,0)=c4p0
c
      cbb4bbx(0,0)=0.d0
      ctp4tpx(0,0)=0.d0
c
      cc0cvx(0,0)=0.d0
      cc2cvx(0,0)=0.d0
      ct2tvx(0,0)=0.d0
      cb2bvx(0,0)=0.d0
c
      cc0ctx(0,0)=0.d0
      cc2ctx(0,0)=0.d0
      ct2ttx(0,0)=0.d0
      cb2btx(0,0)=0.d0
c
c   S=0 and T=1 channel LECs follow
c
      cc0cx(0,1)=c0c0
      cc2cx(0,1)=c2c01
      cc4cx(0,1)=c4c01
c
      ct2tx(0,1)=0.d0
      ct4tx(0,1)=0.d0
c
      cb2bx(0,1)=0.d0
      cb4bx(0,1)=0.d0
c
      cq4qx(0,1)=c4q0
      cp4px(0,1)=c4p0
c
      cbb4bbx(0,1)=0.d0
      ctp4tpx(0,1)=0.d0
c
      cc0cvx(0,1)=dble( 2*itz )*c0cv
      cc2cvx(0,1)=dble( 2*itz )*c2c0v
      ct2tvx(0,1)=0.d0
      cb2bvx(0,1)=0.d0
c
      cc0ctx(0,1)=dble(-4+6*iabs( itz ) )*c0ct
      cc2ctx(0,1)=dble(-4+6*iabs( itz ) )*c2c0t
      ct2ttx(0,1)=0.d0
      cb2btx(0,1)=0.d0
c
c   S=1 and T=0 channel LECs follow
c
      cc0cx(1,0)=c0c1
      cc2cx(1,0)=c2c10
      cc4cx(1,0)=c4c10
c
      ct2tx(1,0)=c2t0
      ct4tx(1,0)=c4t0
c
      cb2bx(1,0)=c2b
      cb4bx(1,0)=c4b0
c
      cq4qx(1,0)=c4q1
      cp4px(1,0)=c4p1
c
      cbb4bbx(1,0)=c4bb
      ctp4tpx(1,0)=c4tp0
c
      cc0cvx(1,0)=0.d0
      cc2cvx(1,0)=0.d0
      ct2tvx(1,0)=0.d0
      cb2bvx(1,0)=0.d0
c
      cc0ctx(1,0)=0.d0
      cc2ctx(1,0)=0.d0
      ct2ttx(1,0)=0.d0
      cb2btx(1,0)=0.d0
c
c   S=1 and T=1 channel LECs follow
c
      cc0cx(1,1)=c0c1
      cc2cx(1,1)=c2c11
      cc4cx(1,1)=c4c11
c
      ct2tx(1,1)=c2t1
      ct4tx(1,1)=c4t1
c
      cb2bx(1,1)=c2b
      cb4bx(1,1)=c4b1
c
      cq4qx(1,1)=c4q1
      cp4px(1,1)=c4p1
c
      cbb4bbx(1,1)=c4bb
      ctp4tpx(1,1)=c4tp1
c
      cc0cvx(1,1)=dble( 2*itz )*c0cv
      cc2cvx(1,1)=dble( 2*itz )*c2c1v
      ct2tvx(1,1)=dble( 2*itz )*c2tv
      cb2bvx(1,1)=dble( 2*itz )*c2bv
c
      cc0ctx(1,1)=dble(-4+6*iabs( itz ) )*c0ct
      cc2ctx(1,1)=dble(-4+6*iabs( itz ) )*c2c1t
      ct2ttx(1,1)=dble(-4+6*iabs( itz ) )*c2tt
      cb2btx(1,1)=dble(-4+6*iabs( itz ) )*c2bt
c
      do 10 it=0,1
      do 10 is=0,1
      cc0cx(is,it)=hbarc*cc0cx(is,it)
      cc2cx(is,it)=hbarc*cc2cx(is,it)
      cc4cx(is,it)=hbarc*cc4cx(is,it)
c
      ct2tx(is,it)=hbarc*ct2tx(is,it)
      ct4tx(is,it)=hbarc*ct4tx(is,it)
c
      cb2bx(is,it)=hbarc*cb2bx(is,it)
      cb4bx(is,it)=hbarc*cb4bx(is,it)
c
      cq4qx(is,it)=hbarc*cq4qx(is,it)
      cp4px(is,it)=hbarc*cp4px(is,it)
c
      cbb4bbx(is,it)=hbarc*cbb4bbx(is,it)
      ctp4tpx(is,it)=hbarc*ctp4tpx(is,it)
c
      cc0cvx(is,it)=hbarc*cc0cvx(is,it)
      cc2cvx(is,it)=hbarc*cc2cvx(is,it)
      ct2tvx(is,it)=hbarc*ct2tvx(is,it)
      cb2bvx(is,it)=hbarc*cb2bvx(is,it)
c
      cc0ctx(is,it)=hbarc*cc0ctx(is,it)
      cc2ctx(is,it)=hbarc*cc2ctx(is,it)
      ct2ttx(is,it)=hbarc*ct2ttx(is,it)
      cb2btx(is,it)=hbarc*cb2btx(is,it)
 10   continue
c
      return
 5000 format(a1,9x,e15.9,1x,a5)
 5010 format(/1x,'rcf  =',e15.9,1x,//,
     &        1x,'c0c0 =',e15.9,1x,'c0c1 =',e15.9,/,
     &        1x,'c2c00=',e15.9,1x,'c2c10=',e15.9,
     &        1x,'c2c01=',e15.9,1x,'c2c11=',e15.9,/,
     &        1x,'c4c00=',e15.9,1x,'c4c10=',e15.9,
     &        1x,'c4c01=',e15.9,1x,'c4c11=',e15.9,/,
     &        1x,'c2t0 =',e15.9,1x,'c2t1 =',e15.9,
     &        1x,'c4t0 =',e15.9,1x,'c4t1 =',e15.9,/,
     &        1x,'c2b  =',e15.9,1x,'c4b0 =',e15.9,
     &        1x,'c4b1 =',e15.9,1x,'c4bb =',e15.9,/,
     &        1x,'c4q0 =',e15.9,1x,'c4q1 =',e15.9,
     &        1x,'c4p0 =',e15.9,1x,'c4p1 =',e15.9,/,
     &        1x,'c4tp0=',e15.9,1x,'c4tp1=',e15.9,/,
     &        1x,'c0cv =',e15.9,1x,'c2c0v=',e15.9,1x,'c2c1v=',e15.9,/,
     &        1x,'c2tv =',e15.9,1x,'c2bv =',e15.9,/
     &        1x,'c0ct =',e15.9,1x,'c2c0t=',e15.9,1x,'c2c1t=',e15.9,/,
     &        1x,'c2tt =',e15.9,1x,'c2bt =',e15.9,//,
     &        1x,'ccc1 =',e15.9,1x,'ccc2 =',e15.9,
     &        1x,'ccc3 =',e15.9,1x,'ccc4 =',e15.9,1x,'cb38 =',e15.9,//,
     &        1x,'chan0 =',e15.9,//)
      end
c
c
c
c_______________________________________________________________________
      subroutine chiral_nlo( idlt , rr , vv , vvd )
c_______________________________________________________________________
c  Potential including Deltas at NLO in r-space obtained using 
c  Krebs et al. (arXiv:nucl-th/0703087v1) and Kaiser et. al 
c  (arXiv:nucl-th/9802071v1) papers.
c
c  idlt = 0 or 1 not to or to compute TPE components including Delta's
c
c  rr = inter-nucleon separation in fm units
c
c  rcf = range in fm units for r-space cutoff 
c
c  vv = OPE and TPE (w/o Delta's) potential components in MeV units:
c
c       vv(1) --> central (it vanishes)               , denoted as c
c       vv(2) --> central x tau1.tau2                 , denoted as ct
c       vv(3) --> central x sigma1.sigma2             , denoted as cs
c       vv(4) --> central x sigma1.sigma2 x tau1.tau2 , denoted as cst
c       vv(5) --> tensor                              , denoted as t 
c       vv(6) --> tensor x tau1.tau2                  , denoted as tt 
c       vv(7) --> sigma1.sigma2 x iso-tensor          , denoted as sit
c       vv(8) --> tensor x iso-tensor                 , denoted as tit
c
c  components vv(4), vv(6) are from OPE only, while
c  components vv(2), vv(3), and vv(5) are from TPE only
c
c  vvd = TPE potential components with Delta's in MeV units:
c
c       vvd(i,1) = component from triangle and box with one Delta diagrams
c       vvd(i,2) = component from box with two Delta's diagrams 
c 
c       for i=1, ...6, i.e., c, ct, cs, cst, t, tt components; note that
c       components vvd(i,1) are proportional to h_A^2, while components
c       vvd(i,2) are proportional to h_A^4, h_A being the N to Delta axial
c       coupling constant, the large N_c value for which usually adopted is
c       h_A = 3 * g_A / sqrt(2)
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      save
c
      include "params_new_pot.f"
      parameter( rx=15.d0 )
c
      common/gf/gx(nz),wgx(nz)
      common/lecs/rcf,c0c0,c0c1,c2c00,c2c10,c2c01,c2c11,
     x                          c4c00,c4c10,c4c01,c4c11,
     x                c2t0,c2t1,c4t0,c4t1,c2b,c4b0,c4b1,c4bb,
     x                c4q0,c4q1,c4p0,c4p1,c4tp0,c4tp1,
     x                c0cv,c2c0v,c2c1v,c2tv,c2bv,
     x                c0ct,c2c0t,c2c1t,c2tt,c2bt,
     x                ccc1,ccc2,ccc3,ccc4,cb38,chan0
c
      dimension vv(8), vvd(6,2)
c
      han0=chan0
c
      do 10 i=1,6
      vv(i)=0.d0
c
      vvd(i,1)=0.d0
      vvd(i,2)=0.d0
 10   continue
c
      vv(7)=0.d0
      vv(8)=0.d0
c
      if( rr .eq. 0.d0 ) return
      if( rr .ge. rx ) go to 30
      alb=xf/rr
c
      albr=alb*rr
c
      xxp =amp*rr
c
      xxp2=xxp *xxp
      xxp3=xxp2*xxp
      xxp4=xxp3*xxp
c
      xxd =amd*rr
      xxd2=xxd*xxd
c
      rr1=1.d0/rr
      rr2=1.d0/rr**2
      rr3=1.d0/rr**3
      rr4=1.d0/rr**4
      rr5=1.d0/rr**5
c
      exp1=dexp(-xxp )*rr1
      exp2=dexp(-xxp )*rr2
      exp3=dexp(-xxp )*rr3
c
      bk0=bkl( 0 , 2.d0*xxp )
      bk1=bkl( 1 , 2.d0*xxp )
c
      ft =aft *rr4*(-12.d0*xxp*bk0-(15.d0+4.d0*xxp2 )*bk1 )
      fcs=afcs*rr4*(  3.d0*xxp*bk0+( 3.d0+2.d0*xxp2 )*bk1 )
      fct=afct*rr4*( ( 1.d0+2.d0*gan2*(  2.d0*xxp2+ 5.d0 )
     x                     -     gan4*( 12.d0*xxp2+23.d0 ) )*bk1
     x              +xxp*( 1.d0+10.d0*gan2
     x                         -      gan4*( 4.d0*xxp2+23.d0 ) )*bk0 )
      
c
      fcst=afcst*exp1
      ftt =aftt *exp3*( 3.d0+3.d0*xxp+xxp2 )      
c
      vv(2)=fct
      vv(3)=fcs
      vv(4)=fcst
      vv(5)=ft
      vv(6)=ftt
c
      if( idlt .eq. 1 ) then
c
      cfctr  =0.d0
      cfctb1 =0.d0
      cfctb2 =0.d0
c      
      cfcb2  =0.d0
c
      cftb1  =0.d0
      cftb2  =0.d0
c
      cfcsb1 =0.d0
      cfcsb2 =0.d0
c
      cfcstb2=0.d0
c
      cfttb2 =0.d0
c
      do 20 im=1,nz
      gim=1.d0-gx(im)
c
      xim =-albr*dlog( gim )
      wxim= albr*wgx(im)/gim

      dp4=xim**2+4.d0*xxp2
      sdp4=dsqrt(dp4)
      exp4=dexp(-sdp4)
c
      cpcb21=-4.d0*xxd2+2.d0*(-2.d0*xxp2-xim**2-2.d0*xxd2)**2
     x                                     /(-xim**2-4.d0*xxd2)
      cpcb22=(-2.d0*xxp2-xim**2-2.d0*xxd2)*(-2.d0*xxp2-xim**2+6.d0*xxd2)
c
      cpctb21=       (-24.d0*xxp2-11.d0*xim**2-24.d0*xxd2)
     x         +6.d0*( -2.d0*xxp2-xim**2-2.d0*xxd2)**2
     x                                    /(-xim**2-4.d0*xxd2)
      cpctb22=3.d0*(-2.d0*xxp2-xim**2-2.d0*xxd2)
     x            *( 2.d0*xxp2+xim**2+10.d0*xxd2)
c
      cpctr11=-12.d0*xxp2-5.d0*xim**2-12.d0*xxd2
      cpctr12= -2.d0*xxp2-     xim**2- 2.d0*xxd2
c
      cpctb11=12.d0*xxd2+24.d0*xxp2+11.d0*xim**2
      cpctb12=-2.d0*xxp2-xim**2-2.d0*xxd2
c
      cptb11=       (3.d0+3.d0*sdp4+dp4)
      cptb12=cptb11*(-xim**2-4.d0*xxd2)
c
      cptb21=       (3.d0+3.d0*sdp4+dp4)
      cptb22=cptb21*(12.d0*xxd2+xim**2 )
c
      cpcsb11=(xim**2+4.d0*xxp2)
      cpcsb12=cpcsb11*(-xim**2-4.d0*xxd2)
c
      cpcsb21=(xim**2+4.d0*xxp2)
      cpcsb22=cpcsb21*(12.d0*xxd2+xim**2 )
c
      cpcstb21=(xim**2+4.d0*xxp2)
      cpcstb22=cpcstb21*(4.d0*xxd2-xim**2 )
c
      cpttb21=        (3.d0+3.d0*sdp4+dp4)
      cpttb22=cpttb21*(4.d0*xxd2-xim**2 )
c
      cfcb21=-xim**2*cpcb21*exp4                          /sdp4
      cfcb22= xim   *cpcb22*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cfctb21=-xim**2*cpctb21*exp4                          /sdp4
      cfctb22= xim   *cpctb22*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cfctr11=-         xim**2*cpctr11*exp4                        /sdp4
      cfctr12= 12.d0*xxd*xim   *cpctr12*exp4*datan( 0.5d0*xim/xxd )/sdp4
c
      cfctb11=-    xim**2*cpctb11*   exp4                         /sdp4
      cfctb12=6.d0*xim   *cpctb12**2*exp4
     x                                  *datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cftb11=2.d0*xim**2*cptb11*exp4                          /sdp4
      cftb12=     xim   *cptb12*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cftb21=-6.d0*xim**2*cptb21*exp4                          /sdp4
      cftb22=      xim   *cptb22*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cfcstb21=-2.d0*xim**2*cpcstb21*exp4                          /sdp4
      cfcstb22=     xim   *cpcstb22*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cfttb21=-2.d0*xim**2*cpttb21*exp4                          /sdp4
      cfttb22=      xim   *cpttb22*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cfcsb11=2.d0*xim**2*cpcsb11*exp4                          /sdp4
      cfcsb12=     xim   *cpcsb12*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cfcsb21=-6.d0*xim**2*cpcsb21*exp4                          /sdp4
      cfcsb22=      xim   *cpcsb22*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cfcb2 =cfcb2+wxim*(cfcb21+cfcb22)
c
      cfctr =cfctr +wxim*(cfctr11+cfctr12)
c
      cfctb1=cfctb1+wxim*(cfctb11+cfctb12)
      cfctb2=cfctb2+wxim*(cfctb21+cfctb22)
c
      cftb1 =cftb1 +wxim*(cftb11+cftb12)
      cftb2 =cftb2 +wxim*(cftb21+cftb22)
c
      cfcsb1=cfcsb1+wxim*(cfcsb11+cfcsb12)
      cfcsb2=cfcsb2+wxim*(cfcsb21+cfcsb22)
c
      cfcstb2=cfcstb2+wxim*(cfcstb21+cfcstb22)
      cfttb2 =cfttb2 +wxim*(cfttb21+cfttb22)
 20   continue
c
      fexp=dexp(-2.d0*xxp)
c
      cfcb1  =fexp*(6.d0+12.d0*xxp+10.d0*xxp2+4.d0*xxp3+xxp4)/xxd
      cfcstb1=fexp*(1.d0+xxp)*(3.d0+3.d0*xxp+2.d0*xxp2      )/xxd
      cfttb1 =fexp*(1.d0+xxp)*(3.d0+3.d0*xxp+     xxp2      )/xxd
c
      ctr=  han0       **2/fpp0**4/pi2 
      cb1=( han0*gan0 )**2/fpp0**4/pi2 
      cb2=( han0/fpp0 )**4        /pi2 
c
      fcb1  =-cb1*rr5*cfcb1 /     6.d0
      fcb2  =-cb2*rr5*cfcb2  /(  108.d0*pi)
      fcttr =-ctr*rr5*cfctr  /(  216.d0*pi)
      fctb1 =-cb1*rr5*cfctb1 /(  216.d0*pi)
      fctb2 =-cb2*rr5*cfctb2 /( 1944.d0*pi)
      ftb1  = cb1*rr5*cftb1  /(  144.d0*pi)
      ftb2  = cb2*rr5*cftb2  /( 2592.d0*pi)
      fcsb1 =-cb1*rr5*cfcsb1 /(   72.d0*pi)
      fcsb2 =-cb2*rr5*cfcsb2 /( 1296.d0*pi)
      fcstb1= cb1*rr5*cfcstb1/    54.d0
      fcstb2=-cb2*rr5*cfcstb2/( 7776.d0*pi)
      fttb1 =-cb1*rr5*cfttb1 /    54.d0
      fttb2 = cb2*rr5*cfttb2 /(15552.d0*pi)
c
      vvd(1,1)=fcb1
      vvd(2,1)=fcttr+fctb1
      vvd(3,1)=fcsb1
      vvd(4,1)=fcstb1
      vvd(5,1)=ftb1
      vvd(6,1)=fttb1
c
      vvd(1,2)=fcb2
      vvd(2,2)=fctb2
      vvd(3,2)=fcsb2
      vvd(4,2)=fcstb2
      vvd(5,2)=ftb2
      vvd(6,2)=fttb2
c
      endif
c
 30   ct20=ct2*ap02/( 12.d0*pi )
      ct2c=ct2*apc2/( 12.d0*pi )
c
      ap0r=rr*ampi0
      apcr=rr*ampic
c 
      yc0=ct20*dexp(-ap0r )/rr
      ycc=ct2c*dexp(-apcr )/rr
c
      yt0=yc0*( 1.d0+3.d0/ap0r+3.d0/ap0r**2 )
      ytc=ycc*( 1.d0+3.d0/apcr+3.d0/apcr**2 ) 
c
      vv(4)=( yc0+2.d0*ycc )/3.d0
      vv(6)=( yt0+2.d0*ytc )/3.d0
c
      vv(7)=( yc0-ycc )/3.d0
      vv(8)=( yt0-ytc )/3.d0
c
c------------------ATTENTION-------------------------
c Choice of r-space cutoff
c
       acf=rcf/2.d0
       wcf=1.d0-1.d0/( ( rr/rcf )**6*dexp( ( rr-rcf )/acf )+1.d0 )
c
c with acf = rcf / 2 
c
      do 60 i=1,6
      vv(i)=hbarc*wcf*vv(i)
c
      vvd(i,1)=hbarc*wcf*vvd(i,1)
      vvd(i,2)=hbarc*wcf*vvd(i,2)
 60   continue
c
      vv(7)=hbarc*wcf*vv(7)
      vv(8)=hbarc*wcf*vv(8)
c
      return
      end
c
c
c
c_______________________________________________________________________
      subroutine chiral_n2lo( idlt , rr , vv2, vvd2 )
c_______________________________________________________________________
c  Potential including Deltas at N2LO in r-space obtained using 
c  Krebs et al. (arXiv:nucl-th/0703087v1) and Epelbaum et. al
c  (Eur. Phys. J. A 19 125-137) papers.
c
c  idlt = 0 or 1 not to or to compute TPE components including Delta's
c
c  rr = inter-nucleon separation in fm units
c
c  rcf = range in fm units for r-space cutoff
c
c  vv2 = TPE (w/o Delta's) potential components at n2lo in MeV units:
c
c       vv2(1) --> central                              , denoted as c
c       vv2(2) --> central x tau1.tau2                  , denoted as ct
c       vv2(3) --> central x sigma1.sigma2              , denoted as cs
c       vv2(4) --> central x sigma1.sigma2 x tau1.tau2  , denoted as cst
c       vv2(5) --> tensor                               , denoted as t 
c       vv2(6) --> tensor x tau1.tau2                   , denoted as tt 
c
c  components vv2(4), vv2(4), vv2(6) are from TPE at n2lo with no Delta's
c
c  vvd2 = TPE potential components with Delta's at n2lo in MeV units:
c
c       vvd2(i,1) = component from triangle and box with one Delta diagrams
c       vvd2(i,2) = component from box with two Delta's diagrams 
c 
c       for i=1, ...6, i.e., c, ct, cs, cst, t, tt components; note that
c       h_A = 3 * g_A / sqrt(2)
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      save
c
      include "params_new_pot.f"
c
      common/gf/gx(nz),wgx(nz)
      common/lecs/rcf,c0c0,c0c1,c2c00,c2c10,c2c01,c2c11,
     x                          c4c00,c4c10,c4c01,c4c11,
     x                c2t0,c2t1,c4t0,c4t1,c2b,c4b0,c4b1,c4bb,
     x                c4q0,c4q1,c4p0,c4p1,c4tp0,c4tp1,
     x                c0cv,c2c0v,c2c1v,c2tv,c2bv,
     x                c0ct,c2c0t,c2c1t,c2tt,c2bt,
     x                ccc1,ccc2,ccc3,ccc4,cb38,chan0
c
      dimension vv2(6),vvd2(6,2)
c
      han0=chan0
c
      xxp =amp*rr
c
      xxp2=xxp *xxp
      xxp3=xxp2*xxp
      xxp4=xxp3*xxp
c
      xxd =amd*rr
      xxd2=xxd*xxd
c
      do 10 i=1,6
      vv2(i)=0.d0
c
      vvd2(i,1)=0.d0
      vvd2(i,2)=0.d0
 10   continue
c
      alb=xf/rr
c
      albr=alb*rr
c
      rr1=1.d0/rr
      rr2=1.d0/rr**2
      rr3=1.d0/rr**3
      rr4=1.d0/rr**4
      rr5=1.d0/rr**5
      rr6=1.d0/rr**6
      
c
      cc1=ccc1*hbarc/1000.d0
      cc2=ccc2*hbarc/1000.d0
      cc3=ccc3*hbarc/1000.d0
      cc4=ccc4*hbarc/1000.d0
      b38=cb38*hbarc/1000.d0
c
      cnd=gan0**2/fpp0**4/pi**2
c
      cndc =   3.d0*cnd      /2.d0
      cndcst= (     cnd*cc4 )/3.d0
      cndtt=- (     cnd*cc4 )/3.d0
c
      fexp=dexp(-2.d0*xxp)
c
      fc  =cndc  *rr6*fexp*( 2.d0*cc1*xxp2*(1.d0+xxp)**2
     x                 +cc3*( 6.d0+12.d0*xxp+10.d0*xxp2+4.d0*xxp3+xxp4))
      fcst=cndcst*rr6*fexp*(1.d0+xxp)*(3.d0+3.d0*xxp+2.d0*xxp2)
      ftt =cndtt *rr6*fexp*(1.d0+xxp)*(3.d0+3.d0*xxp+     xxp2)  
c
      vv2(1)=fc
      vv2(4)=fcst
      vv2(6)=ftt
c
      if( idlt .eq. 1 ) then
c
      cfc2tr=0.d0
      cfc2b2=0.d0
c
      cfc2ttr =0.d0
      cfc2tb1 =0.d0
      cfc2tb2 =0.d0
c      
      cfc2sb1 =0.d0
      cfc2sb2 =0.d0
c
      cfc2sttr =0.d0
      cfc2stb2 =0.d0
c
      cf2tb1=0.d0
      cf2tb2=0.d0
c
      cf2tttr =0.d0
      cf2ttb2 =0.d0
c
      do im=1,nz
      gim=1.d0-gx(im)
c
      xim =-albr*dlog( gim )
      wxim= albr*wgx(im)/gim

      dp4=xim**2+4.d0*xxp2
      sdp4=dsqrt(dp4)
      exp4=dexp(-sdp4)
c
      cpc2tr11=( -24.d0*cc1*xxp2+cc2*(5.d0*xim**2+12.d0*xxp2+12.d0*xxd2)
     x                     -6.d0*cc3*(xim**2+2.d0*xxp2) )
      cpc2tr12=(xim**2+2.d0*xxp2+2.d0*xxd2)*(4.d0*cc1*xxp2-2.d0*cc2*xxd2
     x                                         +cc3*(xim**2+2.d0*xxp2) )
c
      cpc2b21=( 6.d0*(xim**2+2.d0*xxp2+2.d0*xxd2)**2/(xim**2+4.d0*xxd2)
     x          +11.d0*xim**2+24.d0*xxp2+12.d0*xxd2 )
      cpc2b22=(xim**2+2.d0*xxp2+10.d0*xxd2)*(xim**2+2.d0*xxp2+2.d0*xxd2)
c
      cpc2ttr11=( 5.d0*xim**2+12.d0*xxp2+12.d0*xxd2 )
      cpc2ttr12=(      xim**2+ 2.d0*xxp2+ 2.d0*xxd2 )
c
      cpc2tb11=( 11.d0*xim**2+24.d0*xxp2+12.d0*xxd2 )
      cpc2tb12=( xim**2+2.d0*xxp2+2.d0*xxd2 )**2
c
      cpc2tb21=cpc2b21
      cpc2tb22=cpc2b22
c
      cpc2sb11=(xim**2+4.d0*xxp2)
      cpc2sb12=(xim**2+4.d0*xxp2)*(xim**2+4.d0*xxd2)
c
      cpc2sb21=(xim**2+4.d0*xxp2)
      cpc2sb22=(xim**2+4.d0*xxp2)*(xim**2+12.d0*xxd2)
c      
      cpc2sttr11=(xim**2+4.d0*xxp2)
      cpc2sttr12=(xim**2+4.d0*xxp2)*(xim**2+4.d0*xxd2)
c
      cpc2stb21=cpc2sb21
      cpc2stb22=cpc2sb22
c
      cp2tb11=(3.d0+3.d0*sdp4+dp4)     
      cp2tb12=cp2tb11*(xim**2+4.d0*xxd2)
c
      cp2tb21=(3.d0+3.d0*sdp4+dp4) 
      cp2tb22=cp2tb21*(xim**2+12.d0*xxd2)
c
      cp2tttr11=(3.d0+3.d0*sdp4+dp4)
      cp2tttr12=cp2tttr11*(xim**2+4.d0*xxd2)
c
      cp2ttb21=cp2tb21
      cp2ttb22=cp2tb22
c
c      
      cfc2tr11=     xim**2*cpc2tr11*exp4                         /sdp4
      cfc2tr12=6.d0*xim   *cpc2tr12*exp4*datan(0.5d0*xim/xxd)
     x                                                      /(xxd*sdp4)
c
      cfc2b21=      xim**2*cpc2b21*exp4                         /sdp4
      cfc2b22=-3.d0*xim   *cpc2b22*exp4*datan(0.5d0*xim/xxd)
     x                                                     /(xxd*sdp4)
c
      cfc2ttr11=           xim**2*cpc2ttr11*exp4                  /sdp4
      cfc2ttr12=-12.d0*xxd*xim*cpc2ttr12*exp4*datan(0.5d0*xim/xxd)/sdp4
c
      cfc2tb11=-    xim**2*cpc2tb11*exp4                         /sdp4
      cfc2tb12=6.d0*xim   *cpc2tb12*exp4*datan(0.5d0*xim/xxd)
     x                                                      /(xxd*sdp4)
c
      cfc2tb21=      xim**2*cpc2tb21*exp4                         /sdp4
      cfc2tb22=-3.d0*xim   *cpc2tb22*exp4*datan(0.5d0*xim/xxd)
     x                                                       /(xxd*sdp4)
c
      cfc2sb11= 2.d0*xim**2*cpc2sb11*exp4                         /sdp4
      cfc2sb12=-    xim   *cpc2sb12*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cfc2sb21=-6.d0*xim**2*cpc2sb21*exp4                          /sdp4
      cfc2sb22=     xim   *cpc2sb22*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cfc2sttr11=2.d0*xim**2*cpc2sttr11*exp4                       /sdp4
      cfc2sttr12=-    xim   *cpc2sttr12*exp4*datan(0.5d0*xim/xxd)
     x                                                       /(xxd*sdp4)
c
      cfc2stb21=-6.d0*xim**2*cpc2stb21*exp4                        /sdp4
      cfc2stb22=      xim   *cpc2stb22*exp4*datan(0.5d0*xim/xxd)
     x                                                       /(xxd*sdp4)
c
      cf2tb11=2.d0*xim**2*cp2tb11*exp4                            /sdp4
      cf2tb12=-    xim   *cp2tb12*exp4*datan(0.5d0*xim/xxd)  /(xxd*sdp4)
c
      cf2tb21=-6.d0*xim**2*cp2tb21*exp4                          /sdp4
      cf2tb22=      xim   *cp2tb22*exp4*datan(0.5d0*xim/xxd)/(xxd*sdp4)
c
      cf2tttr11=2.d0*xim**2*cp2tttr11*exp4                        /sdp4
      cf2tttr12=-    xim   *cp2tttr12*exp4*datan(0.5d0*xim/xxd)
     x                                                       /(xxd*sdp4)
c
      cf2ttb21=-6.d0*xim**2*cp2ttb21*exp4                        /sdp4
      cf2ttb22=      xim   *cp2ttb22*exp4*datan(0.5d0*xim/xxd)
     x                                                      /(xxd*sdp4)
c     
      cfc2tr=cfc2tr+wxim*(cfc2tr11+cfc2tr12)
      cfc2b2=cfc2b2+wxim*(cfc2b21+cfc2b22)
c
      cfc2ttr=cfc2ttr+wxim*(cfc2ttr11+cfc2ttr12)
      cfc2tb1=cfc2tb1+wxim*(cfc2tb11+cfc2tb12)
      cfc2tb2=cfc2tb2+wxim*(cfc2tb21+cfc2tb22)
c
      cfc2sb1=cfc2sb1+wxim*(cfc2sb11+cfc2sb12)
      cfc2sb2=cfc2sb2+wxim*(cfc2sb21+cfc2sb22)
c 
      cfc2sttr=cfc2sttr+wxim*(cfc2sttr11+cfc2sttr12)
      cfc2stb2=cfc2stb2+wxim*(cfc2stb21+cfc2stb22)
c
      cf2tb1=cf2tb1+wxim*(cf2tb11+cf2tb12)
      cf2tb2=cf2tb2+wxim*(cf2tb21+cf2tb22)
c
      cf2tttr=cf2tttr+wxim*(cf2tttr11+cf2tttr12)
      cf2ttb2=cf2ttb2+wxim*(cf2ttb21+cf2ttb22)
      enddo
c
      ctr2 =han0**2/fpp0**4/pi**3
      ctr38=han0   /fpp0**4/pi**3
      cb12=han0*gan0**2/fpp0**4/pi**3
      cb22=han0**3/fpp0**4/pi**3

      fc2tr=          ctr2*rr6*xxd*cfc2tr/18.d0
      fc2b2=-2.d0*b38*cb22*rr6*xxd*cfc2b2/81.d0
c
      fc2ttr=-    ctr38*b38*rr6*xxd*cfc2ttr/54.d0
      fc2tb1=-    cb12 *b38*rr6*xxd*cfc2tb1/54.d0
      fc2tb2=-2.d0*cb22*b38*rr6*xxd*cfc2tb2/486.d0
c
      fc2sb1=-    cb12 *b38*rr6*xxd*cfc2sb1/18.d0
      fc2sb2=-    cb22 *b38*rr6*xxd*cfc2sb2/162.d0
c
      fc2sttr=-ctr2*cc4*rr6*xxd*cfc2sttr/108.d0
      fc2stb2=-    cb22 *b38*rr6*xxd*cfc2stb2/972.d0
c
      f2tb1=cb12 *b38*rr6*xxd*cf2tb1/36.d0
      f2tb2=cb22 *b38*rr6*xxd*cf2tb2/324.d0
c
      f2tttr=ctr2*cc4*rr6*xxd*cf2tttr/216.d0
      f2ttb2=cb22 *b38*rr6*xxd*cf2ttb2/1944.d0
c     
      vvd2(1,1)=fc2tr
      vvd2(2,1)=fc2ttr+fc2tb1
      vvd2(3,1)=fc2sb1
      vvd2(4,1)=fc2sttr
      vvd2(5,1)=f2tb1
      vvd2(6,1)=f2tttr
c
      vvd2(1,2)=fc2b2
      vvd2(2,2)=fc2tb2
      vvd2(3,2)=fc2sb2
      vvd2(4,2)=fc2stb2
      vvd2(5,2)=f2tb2
      vvd2(6,2)=f2ttb2
c 
      endif
c
c------------------ATTENTION-------------------------
c Choice of r-space cutoff
c
      acf = rcf / 2.d0
      wcf=1.d0-1.d0/( (rr/rcf)**6*dexp( ( rr-rcf )/acf )+1.d0 )
c 
      do 30 i=1,6
      vv2(i)=hbarc*wcf*vv2(i)
c
      vvd2(i,1)=hbarc*wcf*vvd2(i,1)
      vvd2(i,2)=hbarc*wcf*vvd2(i,2)
 30   continue
c
      return
      end
c
c
c_______________________________________________________________________
      subroutine fnt_cntc( n , ilb , rr ,
     x                     fc0  , fc2   , ft2 , fb2  ,
     x                     fc4  , ft4   , fb4 , fbb4 , fq4 ,
     x                     fp4  , fp4d  ,
     x                     ftp4 , ftp4d )
c_______________________________________________________________________
c
c Contact radial functions at LO, NLO, N3LO:
c
c n = exponent for r-space cutoff
c
c rcf = range in fm units for r-space cutoff
c
c rr = inter-nucleon separation in fm units
c
c fc0   = v^c(r) at LO in fm^{-3} units
c fc2   = v^c(r) at NLO in fm^{-5} units
c ft2   = v^t(r) at NLO in fm^{-5} units
c fb2   = v^b(r) at NLO in fm^{-5} units
c fc4   = v^c(r) at N3LO in fm^{-7} units
c ft4   = v^t(r) at N3LO in fm^{-7} units
c fb4   = v^b(r) at N3LO in fm^{-7} units
c fbb4  = v^bb(r) at N3LO in fm^{-7} units
c fq4   = v^{q}(r) an N3LO in fm^{-7} units
c fp4   = v^{p}(r) an N3LO in fm^{-5} units
c fp4d  = first derivative of v^{p}(r) an N3LO in fm^{-6} units
c ftp4  = v^{tp}(r) an N3LO in fm^{-5} units
c ftp4d = first derivative of v^{tp}(r) an N3LO in fm^{-6} units
c
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      include "params_new_pot.f"
c
      if ( ilb .eq. 1 .or. ilb .eq. 7 .or. ilb .eq. 8 ) then
      rcf1=0.6d0
c
      elseif ( ilb .eq.  2 .or. ilb .eq. 4  .or.
     x         ilb .eq.  5 .or. ilb .eq. 9  .or. 
     x         ilb .eq. 11 .or. ilb .eq. 12 .or.
     x         ilb .eq. 13 ) then
      rcf1=0.7d0
c
      elseif ( ilb .eq. 3 .or. ilb .eq. 6 .or. ilb .eq. 10 ) then
      rcf1=0.8d0
      endif
c
      rcfn=rcf1**n
c
      n1=n-1
      n2=n-2
      n3=n-3
      n4=n-4
c
      rrn0=rr**n
      rrn1=rr**n1
      rrn2=rr**n2
      rrn3=rr**n3
      rrn4=rr**n4
c
      rr2=rr*rr
      rr3=rr*rr2
      rr4=rr*rr3
c
      rcff=dble( n )/rcfn
c
      ann=0.25d0*n/( pi*dgamma( 3.d0/dble( n ) )*rcf1**3 )
      wcf=ann*dexp(-rrn0/rcfn )
c
      wcfd1=-rcff*wcf*rrn1                     
c
      wcfd2=-rcff*(            wcfd1*rrn1
     x             +dble( n1 )*wcf  *rrn2 )
c
      wcfd3=-rcff*(                    wcfd2*rrn1
     x             +2.d0*dble( n1    )*wcfd1*rrn2
     x             +     dble( n1*n2 )*wcf  *rrn3 )
c
      wcfd4=-rcff*(                       wcfd3*rrn1
     x             +3.d0*dble( n1       )*wcfd2*rrn2
     x             +3.d0*dble( n1*n2    )*wcfd1*rrn3
     x             +     dble( n1*n2*n3 )*wcf  *rrn4 )
c
      fc0=wcf
c
      fc2=-( wcfd2+2.d0*wcfd1/rr )
      ft2=-( wcfd2-     wcfd1/rr )
      fb2=-             wcfd1/rr
c
      fc4=  ( wcfd4+4.d0*wcfd3/rr                               )
      ft4=  ( wcfd4+     wcfd3/rr-6.d0*wcfd2/rr2+6.d0*wcfd1/rr3 )
      fb4=  (            wcfd3/rr+2.d0*wcfd2/rr2-2.d0*wcfd1/rr3 )
      fbb4=-(                          wcfd2/rr2-     wcfd1/rr3 )
      fq4 =-(                          wcfd2/rr2-     wcfd1/rr3 )
c
      fp4 =fc2
      ftp4=ft2
c
      fp4d =-wcfd3-2.d0*wcfd2/rr+2.d0*wcfd1/rr2
      ftp4d=-wcfd3+     wcfd2/rr-     wcfd1/rr2
c
      return
      end
c
c
c
c_______________________________________________________________________
      subroutine em_pot( lemp , rr , vvem )
c_______________________________________________________________________
c
c lemp: 0=Coulomb w/ff (pp) only in MeV units
c       1=full electromagnetic potential in MeV units
c
c order of operators in vvem(l)
c l:     1=Coulomb w/ff (pp)   2=DF    (pp)          3=OO      (pp)
c        4=VP    (pp)                                5=Coulomb (np)
c        6=s1.s2 (pp)          7=s1.s2 (nn)          8=s1.s2   (np)
c        9=S12   (pp)         10=S12   (nn)         11=S12     (np)
c       12=L.S   (pp)         13=L.S   (nn)         14=L.S     (np)
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      real*8 me,mp,mn,mr,mup,mun
      save
c 
      include "params_new_pot.f"
      parameter ( alpha=1.d0/137.03599d0 )
c
      dimension vvem(14)
c
      data me / 0.510999d0 / , mp / 938.2720d0 / , mup / 2.79285d0 / ,
     x                         mn / 939.5653d0 / , mun /-1.91304d0 / ,
     x     bb / 4.27d0 / , beta / 0.0189d0 /
c
      br=bb*rr
      fcoul=1.d0-( 1.d0+11.d0*br   /16.d0
     x                 + 3.d0*br**2/16.d0
     x                 +      br**3/48.d0 )*dexp(-br )
      vvem(1)=alpha*hbarc*fcoul/rr
c
      if( lemp .ge. 1 ) then
c
      mr=mp*mn/( mp+mn )
      ha=alpha*hbarc**3
c
      ftensor=1.d0-( 1.d0+br+br**2/  2.d0
     x                      +br**3/  6.d0
     x                      +br**4/ 24.d0
     x                      +br**5/144.d0 )*dexp(-br )
      fspinor=1.d0-( 1.d0+br+     br**2/ 2.d0
     x                      +7.d0*br**3/48.d0
     x                      +     br**4/48.d0 )*dexp(-br )
      fdelta=bb**3*( 1.d0+br+br**2/3.d0 )*dexp(-br )/16.d0
      fnp   =bb**2*( 15.d0*br   +15.d0*br**2
     x              + 6.d0*br**3+      br**4 )*dexp(-br )/384.d0
      fivp=ffvp( 2.d0*me*rr/hbarc )
c-------------------------------------------
c  central pp cpmponents follow
c
      vvem(2)=-ha*fdelta/( 4.d0*mp**2 )
      vvem(3)=-vvem(1)**2/mp
      vvem(4)=2.d0*alpha*vvem(1)*fivp/( 3.d0*pi )
c-------------------------------------------
c  central np component follows
c
      vvem(5)=hbarc*alpha*beta*fnp/rr
c-------------------------------------------
c  s1.s2 components for pp, nn, and np follow
c
      vvem(6)=-ha*mup*mup*fdelta/( 6.d0*mp*mp )
      vvem(7)=-ha*mun*mun*fdelta/( 6.d0*mn*mn )
      vvem(8)=-ha*mup*mun*fdelta/( 6.d0*mn*mp )
c-------------------------------------------
c  tensor components for pp, nn, and np follow
c
      vvem( 9)=-ha*mup*mup*ftensor/( 4.d0*mp*mp*rr**3 )
      vvem(10)=-ha*mun*mun*ftensor/( 4.d0*mn*mn*rr**3 )
      vvem(11)=-ha*mup*mun*ftensor/( 4.d0*mp*mn*rr**3 )
c-------------------------------------------
c  spin-orit components for pp, nn, and np follow
c
      vvem(12)=-ha*( 4.d0*mup-1.d0 )*fspinor/( 2.d0*mp*mp*rr**3 )
      vvem(13)=0.d0
      vvem(14)=-ha*       mun       *fspinor/( 2.d0*mn*mr*rr**3 )
c
      else
c
      do 10 i=2,14
      vvem(i)=0.d0
 10   continue
c
      endif
c
      return
      end
c
c
c
c_______________________________________________________________________
      real*8 function ffvp( zz )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      parameter ( nx=201 , zinf=0.01d0 )
c
      common/gp/gg(nx),wgg(nx)
c
      dimension xx(nx),wx(nx)
c
      data gamma / 0.577215664901533d0 / ,
     x        pi / 3.141592653589793d0 /
c
      ffvp=0.d0
c
      if( zz .le. zinf ) then
      ffvp=-gamma-5.d0/6.d0-dlog( 0.50d0*zz )+3.d0*pi*zz/8.d0
c
      else
c
      alf=1.0d0/zz
      bet=0.5d0*pi*alf
c
      do 10 i=1,nx
      xxi=0.5d0*pi*gg(i)
c
      xx(i)=1.d0+alf*dtan( xxi )
      wx(i)=     bet*wgg(i)/dcos( xxi )**2
 10   continue
c
      do 20 i=1,nx
      x1=xx(i)
      x2=xx(i)*xx(i)
c
      fx=dexp(-zz*( x1-1.d0 ) )*dsqrt( x2-1.d0 )*( 1.d0+0.5d0/x2 )/x2
c
      ffvp=ffvp+wx(i)*fx
 20   continue
c
      ffvp=dexp(-zz )*ffvp
      endif
c
      return
      end
c_______________________________________________________________________
      subroutine axb( aa , bb , cc )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      dimension aa(2,2),bb(2,2),cc(2,2)
c
      do 10 j=1,2
      do 10 i=1,2
      cc(i,j)=0.d0
c
      do 10 l=1,2
      cc(i,j)=cc(i,j)+aa(i,l)*bb(l,j)
 10   continue
c
      return
      end
c
c
c
c_______________________________________________________________________
      subroutine caxcb( ca , cb , cc )
c_______________________________________________________________________
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      save
c
      dimension ca(2,2),cb(2,2),cc(2,2)
c
      do 10 j=1,2
      do 10 i=1,2
      cc(i,j)=dcmplx( 0.d0 , 0.d0 )
c
      do 10 l=1,2
      cc(i,j)=cc(i,j)+ca(i,l)*cb(l,j)
 10   continue
c
      return
      end
c
c
c
c_______________________________________________________________________
      subroutine aiv( aa , av )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      parameter( eps=1.d-32 )
c
      dimension aa(2,2),av(2,2)
c
      dd=aa(1,1)*aa(2,2)-aa(1,2)*aa(2,1)
      if( dabs( dd ) .lt. eps ) go to 9000
c
      av(1,1)= aa(2,2)/dd
      av(2,1)=-aa(2,1)/dd
      av(1,2)=-aa(1,2)/dd
      av(2,2)= aa(1,1)/dd
c
      return
 9000 write(6,5000)
      return
c
 5000 format(10x,//,'determinant smaller than eps',//)
      end
c
c
c
c-----------------------------------------------------------------------
      real*8 function xcb( ll , ml , ms , jj , mj )
c-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      save
c
      xcb=0.d0
c
      if( ll .lt. 0 .or. jj .lt. 0 ) return
      if( ml+ms .ne. mj ) return
      if( jj .lt. iabs( ll-1 ) .or. jj .gt. ( ll+1 ) ) return
c
      if( iabs( ml ) .gt. ll ) return
      if( iabs( ms ) .gt.  1 ) return
      if( iabs( mj ) .gt. jj ) return
c
      if( jj .eq. ll-1 ) then
c
      if( ms .eq. -1 ) then
      xcb=dsqrt( dble( ( ll+mj+1 )*( ll+mj ) )
     x          /dble( 2*ll*( 2*ll+1 ) ) )
      return
c
      elseif( ms .eq. 0 ) then
      xcb=-dsqrt( dble( ( ll-mj )*( ll+mj ) )
     x           /dble( ll*( 2*ll+1 ) ) )
      return
c
      elseif( ms .eq. 1 ) then
      xcb=dsqrt( dble( ( ll-mj+1 )*( ll-mj ) )
     x          /dble( 2*ll*( 2*ll+1 ) ) )
      return
c
      endif
c
      elseif( jj .eq. ll ) then
c
      if( ms .eq. -1 ) then
      xcb=dsqrt( dble( ( ll+mj+1 )*( ll-mj ) )
     x          /dble( 2*ll*( ll+1 ) ) )
      return
c
      elseif( ms .eq. 0 ) then
      xcb=mj/dsqrt( dble( ll*( ll+1 ) ) )
      return
c
      elseif( ms .eq. 1 ) then
      xcb=-dsqrt( dble( ( ll-mj+1 )*( ll+mj ) )
     x           /dble( 2*ll*( ll+1 ) ) )
      return
c
      endif
c
      elseif( jj .eq. ll+1 ) then
c
      if( ms .eq. -1 ) then
      xcb=dsqrt( dble( ( ll-mj+1 )*( ll-mj ) )
     x          /dble( ( 2*ll+2 )*( 2*ll+1 ) ) )
      return
c
      elseif( ms .eq. 0 ) then
      xcb=dsqrt( dble( ( ll-mj+1 )*( ll+mj+1 ) )
     x          /dble( ( 2*ll+1 )*( ll+1 ) ) )
      return
c
      elseif( ms .eq. 1 ) then
      xcb=dsqrt( dble( ( ll+mj+1 )*( ll+mj ) )
     x          /dble( ( 2*ll+2 )*( 2*ll+1 ) ) )
      endif
c
      endif
c
      return
      end
c
c
c
c_______________________________________________________________________________
      real*8 function bnl( l , x )
c_______________________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c     
      blm=-dcos(x)/x
      bl0=-dcos(x)/x**2-dsin(x)/x
c     
      if( l .eq. 0 ) bnl=blm
      if( l .eq. 1 ) bnl=bl0
      if( l .eq. 0 .or. l .eq. 1 ) return
c     
      do 10 k=1,l-1
      blp=dble( 2*k+1 )*bl0/x-blm
      blm=bl0
      bl0=blp
10    continue
c
      bnl=bl0
c     
      return
      end
c
c
c
c________________________________________________________________________
      real*8 function bjl( l , x )
c________________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      data acc / 50.d0 / ,  fac / 1.d+10 /
c
      bjl=dsin(x)/x**( l+1 )-dble(l)*dcos(x)/x
c
      if( l .eq. 0 .or. l .eq. 1 ) return
c
      if( x .gt. dble( l ) ) then
c
      blm=dsin(x)/x
      bl0=dsin(x)/x**2-dcos(x)/x
c
      do 10 k=1,l-1
      blp=dble( 2*k+1 )*bl0/x-blm
      blm=bl0
      bl0=blp
10    continue
c
      bjl=bl0
c
      endif
c
      if( x .le. dble(l) ) then
c
      lsup=2*idint( 0.5d0*( dble(l)+dsqrt( acc*dble(l) ) ) )
      iflag=1
      sume=0.d0
      sumo=0.d0
c
      blp=0.d0
      bl0=1.d0
      bjl=0.d0
c
      do 50 k=lsup,1,-1
      blm=dble( 2*k+1 )*bl0/x-blp
      blp=bl0
      bl0=blm
c
      if( dabs( bl0 ) .gt. fac ) then
c
      blp=blp/fac
      bl0=bl0/fac
      bjl=bjl/fac
      sume=sume/fac
      sumo=sumo/fac
c
      endif
c
      go to(20,30) iflag
c
20    sgnk=dsin( 0.5d0*dble( k-1 )*dacos(-1.d0 ) )
      sumo=sumo+dble( 2*k-1 )*sgnk*bl0
c
      go to 40
c
30    sgnk=dcos( 0.5d0*dble( k-1 )*dacos(-1.d0 ) )
      sume=sume+dble( 2*k-1 )*sgnk*bl0
c
40    iflag=3-iflag
c
      if( k .eq. l ) bjl=blp
50    continue
c
      if( dcos(x) .ne. 0.d0 ) an=sume/dcos(x)
      if( dcos(x) .eq. 0.d0 ) an=sumo/dsin(x)
      bjl=bjl/an
c
      endif
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
      real*8 function fmf( za , ra , el , ml )
c-----------------------------------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      real*8 ml
      save
c
c   This subroutine generates the Fermi function for a nucleus
c   of atomic number ZA and radius RA (in fm); EL and
c   ML are the energy and mass of the lepton (in MeV).
c
      data hc / 197.32696d0 / , alpha / 137.03599d0 / ,
     x     pi / 3.141592653589793d0 /
c
      xx=dsqrt( 1.d0-( za/alpha )**2 )
      yy=za/dsqrt( 1.d0-( ml/el )**2 )/alpha
c
      cc=dcmplx( 1.0d0-xx ,-yy )
      cz=cc+1.0d0
c
      cng=cdlog( pi*cc )-cdlog( cdsin( pi*cc ) )-cgmln( cz )
      cng=cdexp( cng )
c
      cz=dcmplx( 1.d0+2.d0*xx , 0.d0 )
c
      cdg=cgmln( cz )
      cdg=cdexp( cdg )
c
      an2=dreal( cng*dconjg( cng ) )
      ad2=dreal( cdg*dconjg( cdg ) )
      ar=an2/ad2
c
      pl=dsqrt( el**2-ml**2 )/hc
c
      fmf=2.d0*( 1.d0+xx )*( 2.d0*pl*ra )**( 2.d0*xx-2.d0 )
     x                    *dexp( pi*yy )*ar
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
      complex*16 function cgmln( cz )
c-----------------------------------------------------------------------
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      save
c
      dimension pf(6)
c
      data pf / 76.18009173d0  ,-86.50532033d0   , 24.01409822d0 ,
     x          -1.231739516d0 ,  0.120858003d-2 ,  -0.536382d-5 /
      data sp / 2.50662827465d0 /
c
      cc=cz-1.0d0
      cp=cc+5.5d0
c
      cp=( cc+0.5d0 )*cdlog( cp )-cp
c
      cs=dcmplx( 1.d0 , 0.d0 )
c
      do 10 j=1,6
      cc=cc+1.0d0
c
      cs=cs+pf(j)/cc
10    continue
c
      cgmln=cp+cdlog( sp*cs )
c
      return
      end
c
c
c
c_______________________________________________________________________
      subroutine spline( n , x , y , u , yp1 , ypn , y2 )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      dimension x(n),y(n),y2(n),u(n)
c
      if( yp1 .gt. 0.99e30 ) then
      y2(1)=0.
      u (1)=0.
c
      else
c
      y2(1)=-0.5
      u (1)=( 3./( x(2)-x(1) ) )*( ( y(2)-y(1) )/
     x                             ( x(2)-x(1) )-yp1 )
      endif
c
      do 10 i=2,n-1
      sig=( x(i)-x(i-1) )/( x(i+1)-x(i-1) )
      p  =sig*y2(i-1)+2.
c
      y2(i)=( sig-1. )/p
      u (i)=( 6.*( ( y(i+1)-y(i) )/( x(i+1)-x(i) )
     x            -( y(i)-y(i-1) )/( x(i)-x(i-1) ) )
     x       /( x(i+1)-x(i-1) )-sig*u(i-1) )/p
10    continue
c
      if( ypn .gt. 0.99e30 ) then
      qn=0.
      un=0.
c
      else
c
      qn=0.5
      un=( 3./( x(n)-x(n-1) ) )*( ypn-( y(n)-y(n-1) )
     x                               /( x(n)-x(n-1) ) )
      endif
c
      y2(n)=( un-qn*u(n-1) )/( qn*y2(n-1)+1. )
c
      do 20 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
20    continue
c
      return
      end
c
c
c
c_______________________________________________________________________
      subroutine splint( n , xa , ya , y2a , x , y )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      dimension xa(n),ya(n),y2a(n)
c
      klo=1
      khi=n
c
10    if( khi-klo .gt. 1 ) then
      k=( khi+klo )/2
c
      if( xa(k) .gt. x ) then
      khi=k
c
      else
c
      klo=k
      endif
c
      goto 10
      endif
c
      h=xa(khi)-xa(klo)
c
      if( h .eq. 0. ) write(7,5000) 
c
      a=( xa(khi)-x )/h
      b=(-xa(klo)+x )/h
c
      y=a*ya(klo)+b*ya(khi)+( a*( a*a-1. )*y2a(klo)
     x                       +b*( b*b-1. )*y2a(khi) )*h*h/6.
c
      return
5000  format('  bad xa input')
      end
c
c
c
c_________________________________________________________________________
      real*8 function sgm( l , eta )
c_________________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      parameter( eps=1.d-15 )
c
      data gam / 0.577215664901533d0 /
c
      n=0
c
      sgm0=-eta*gam
c
10    etn=eta/dble( 1+n )
      sgm=sgm0+( etn-datan( etn ) )
c
      if( dabs( sgm/sgm0-1.d0 ) .lt. eps ) go to 20
      n=n+1
      sgm0=sgm
      go to 10
c
20    if( l .eq. 0 ) return
c
      do 30 m=0,l-1
      sgm=sgm+datan( eta/dble( 1+m ) )
30    continue
c
      return
      end
c
c
c
c_________________________________________________________________________
      real*8 function sgm0( l , eta )
c_________________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      sgm0=0.d0
      if( l .eq. 0 ) return
c
      do 30 m=1,l
      sgm0=sgm0+datan( eta/dble( m ) )
30    continue
c
      return
      end
c
c
c
c________________________________________________________________________
      real*8 function bkl( l , x )
c________________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c     
      bm=bk0(x)
      b0=bk1(x)
c     
      if( l .le. 1 ) then
      bkl=dble( 1-l )*bm+dble( l )*b0
      else 
c     
      do 10 j=1,l-1
      bp=bm+dble( 2*j )*b0/x
c     
      bm=b0
      b0=bp
10    continue
c
      bkl=b0
      endif
c
      return
      end
c
c
c
c________________________________________________________________________
      real*8 function erfc( x )
c________________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      z=dabs( x )
c
      t=1.d0/( 1.d0+0.5d0*z )
c
      erfc=t*dexp(-z*z-1.26551223d0+t*( 1.00002368d0
     x                             +t*( 0.37409196d0
     x                             +t*( 0.09678418d0
     x                             +t*(-0.18628806d0
     x                             +t*( 0.27886807d0
     x                             +t*(-1.13520398d0
     x                             +t*( 1.48851587d0
     x                             +t*(-0.82215223d0
     x                             +t*  0.17087277d0 ) ) ) ) ) ) ) ) )
c
      if( x .lt. 0.d0 ) erfc=2.d0-erfc
c
      return
      end
c
c
c
c________________________________________________________________________
      real*8 function bk0( x )
c________________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      data p1 /-0.57721566d0 / , p2 / 0.42278420d0 / ,
     x     p3 / 0.23069756d0 / , p4 / 0.3488590d-1 / ,
     x     p5 / 0.262698d-2  / , p6 / 0.10750d-3   / ,
     x     p7 / 0.74d-5 /
c
      data q1 / 1.25331414d0 / , q2 /-0.7832358d-1 / ,
     x     q3 / 0.2189568d-1 / , q4 /-0.1062446d-1 / ,
     x     q5 / 0.587872d-2  / , q6 /-0.251540d-2  / ,
     x     q7 / 0.53208d-3   /
c
      if( x .le. 2.d0 ) then
      y=x*x/4.d0
c
      bk0= ( p1+y*( p2+y*( p3+y*( p4+y*( p5+y*( p6+y*p7 ) ) ) ) ) )
     x     -dlog( x/2.d0 )*bi0( x )
      else
      y=2.d0/x
c
      bk0= ( q1+y*( q2+y*( q3+y*( q4+y*( q5+y*( q6+y*q7 ) ) ) ) ) )
     x    *dexp(-x )/dsqrt( x )
      endif
c
      return
      end
c
c
c
c________________________________________________________________________
      real*8 function bk1( x )
c________________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      data p1 / 1.d0         / , p2 / 0.15443144d0 / ,
     x     p3 /-0.67278579d0 / , p4 /-0.18156897d0 / ,
     x     p5 /-0.1919402d-1 / , p6 /-0.110404d-2  / ,
     x     p7 /-0.4686d-4    /
      data q1 / 1.25331414d0 / , q2 / 0.23498619d0 / ,
     x     q3 /-0.3655620d-1 / , q4 / 0.1504268d-1 / ,
     x     q5 /-0.780353d-2  / , q6 / 0.325614d-2  / ,
     x     q7 /-0.68245d-3   /
c
      if( x .le. 2.d0 ) then
      y=x*x/4.d0
c
      bk1= ( p1+y*( p2+y*( p3+y*( p4+y*( p5+y*( p6+y*p7 ) ) ) ) ) )/x
     x     +dlog( x/2.d0 )*bi1( x )
      else
      y=2.d0/x
c
      bk1= ( q1+y*( q2+y*( q3+y*( q4+y*( q5+y*( q6+y*q7 ) ) ) ) ) )
     x    *dexp(-x )/dsqrt( x )
      endif
c
      return
      end
c
c
c
c________________________________________________________________________
      real*8 function bi0( x )
c________________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      data p1 / 1.0d0       / , p2 / 3.5156229d0 / ,
     x     p3 / 3.0899424d0 / , p4 / 1.2067492d0 / ,
     x     p5 / 0.2659732d0 / , p6 / 0.360768d-1 / ,
     x     p7 / 0.45813d-2  /
c
      data q1 / 0.39894228d0 / , q2 / 0.1328592d-1 / ,
     x     q3 / 0.225319d-2  / , q4 /-0.157565d-2  / ,
     x     q5 / 0.916281d-2  / , q6 /-0.2057706d-1 / ,
     x     q7 / 0.2635537d-1 / , q8 /-0.1647633d-1 / ,
     x     q9 / 0.392377d-2  /
c
      if( dabs( x ) .lt. 3.75d0 ) then
      y=( x/3.75d0 )**2
c
      bi0= p1+y*( p2+y*( p3+y*( p4+y*( p5+y*( p6+y*p7 ) ) ) ) )
      else
      y=3.75d0/dabs( x )
c
      bi0= ( q1+y*( q2+y*( q3+y*( q4+y*( q5+y*( q6
     x                       +y*( q7+y*( q8+y*q9 ) ) ) ) ) ) ) )
     x    *dexp( dabs( x ) )/dsqrt( dabs( x ) )
      endif
c
      return
      end
c
c
c
c________________________________________________________________________
      real*8 function bi1( x )
c________________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      data p1 / 0.5d0        / , p2 / 0.87890594d0 / ,
     x     p3 / 0.51498869d0 / , p4 / 0.15084934d0 / ,
     x     p5 / 0.2658733d-1 / , p6 / 0.301532d-2  / ,
     x     p7 / 0.32411d-3   /
c
      data q1 / 0.39894228d0 / , q2 /-0.3988024d-1 / ,
     x     q3 /-0.362018d-2  / , q4 / 0.163801d-2  / ,
     x     q5 /-0.1031555d-1 / , q6 / 0.2282967d-1 / ,
     x     q7 /-0.2895312d-1 / , q8 / 0.1787654d-1 / ,
     x     q9 /-0.420059d-2  /
c
      if( dabs( x ) .lt. 3.75d0 ) then
      y=( x/3.75d0 )**2
c
      bi1= x*( p1+y*( p2+y*( p3+y*( p4+y*( p5+y*( p6+y*p7 ) ) ) ) ) )
      else
      y=3.75d0/dabs( x )
c
      bi1= ( q1+y*( q2+y*( q3+y*( q4+y*( q5+y*( q6
     x                       +y*( q7+y*( q8+y*q9 ) ) ) ) ) ) ) )
     x    *dexp( dabs( x ) )/dsqrt( dabs( x ) )
      endif
c
      return
      end
c
c
c
c_______________________________________________________________________
      function bode7( ninit , nfin , h , g )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      dimension g(ninit:nfin)
c
      ncount=0
      sum=0.
      nleft=nfin+1-ninit
c
      if( nleft .lt. 8 ) go to 20
c
      n7=7*int( nleft/7 )+ninit-1
      wp=7.d0*h/17280.d0
c
      do 10 i=ninit,n7,7
      ncount=ncount+1
c
      w0= 751.d0*g(i)
      w1=3577.d0*g(i+1)
      w2=1323.d0*g(i+2)
      w3=2989.d0*g(i+3)
      w4=2989.d0*g(i+4)
      w5=1323.d0*g(i+5)
      w6=3577.d0*g(i+6)
      w7= 751.d0*g(i+7)
c
      sum=sum+wp*( w0+w1+w2+w3+w4+w5+w6+w7 )
10    continue
c
      nleft=nfin-7*ncount-ninit+1
c
20    if( nleft .eq. 0 ) go to 90
      n=nfin-nleft+1
c
      go to(90,80,70,60,50,40,30) nleft
30    wp=h/140.d0
c
      wa= 41.d0*( g(n  )+g(n+6) )
      wb=216.d0*( g(n+1)+g(n+5) )
      wc= 27.d0*( g(n+2)+g(n+4) )
      wd=272.d0*  g(n+3)
      sum=sum+wp*( wa+wb+wc+wd )
c
      go to 90
c
40    wp=5.d0*h/288.d0
c
      wa=19.d0*( g(n  )+g(n+5) )
      wb=75.d0*( g(n+1)+g(n+4) )
      wc=50.d0*( g(n+2)+g(n+3) )
      sum=sum+wp*( wa+wb+wc )
c
      go to 90
c
50    wp=2.d0*h/45.d0
c
      wa= 7.d0*( g(n  )+g(n+4) )
      wb=32.d0*( g(n+1)+g(n+3) )
      wc=12.d0*  g(n+2)
      sum=sum+wp*( wa+wb+wc )
c
      go to 90
c
60    wp=3.d0*h/8.d0
c
      wa=       g(n  )+g(n+3)
      wb=3.d0*( g(n+1)+g(n+2) )
      sum=sum+wp*( wa+wb )
c
      go to 90
c
70    wp=h/3.d0
c
      wa=     g(n  )+g(n+2)
      wb=4.d0*g(n+1)
      sum=sum+wp*( wa+wb )
c
      go to 90
c
80    wp=h/2.d0
c
      wa=g(n)+g(n+1)
      sum=sum+wp*wa
c
90    bode7=sum
c
      return
      end
c
c
c
c_______________________________________________________________________
      subroutine gauss_pot( ng , ainf , asup , xg , wg )
c_______________________________________________________________________
      implicit real*8(a-h,o-z)
      save
c
      dimension xg(ng),wg(ng)
c
      data eps / 1.d-12 / 
c
      apl=.5d0*( asup+ainf )
      amn=.5d0*( asup-ainf )
c
      c=.5d0*dacos(-1.d0)/dble( 2*ng+1 )
c
      do 30 i=1,( ng+1 )/2
      x=dcos( c*dble( 4*i-1 ) )
c
10    p0=1.d0
      p1=0.d0
c
      do 20 j=1,ng
      c1=2.d0-1.d0/dble(j)
      c2=1.d0-1.d0/dble(j)
      p2=p1
      p1=p0
      p0=x*c1*p1-c2*p2
20    continue
c
      dp0=dble( ng )*( x*p0-p1 )/( x**2-1.d0 )
      x0=x
      x=x0-p0/dp0
c
      xchk=x-x0
      if( xchk .lt. 0.d0 ) xchk=-xchk
      if( xchk .gt. eps ) go to 10
c
      xg(i     )=apl-amn*x
      xg(ng+1-i)=apl+amn*x
      wg(i     )=2.d0*amn/( (1.d0-x**2)*dp0**2 )
      wg(ng+1-i)=2.d0*amn/( (1.d0-x**2)*dp0**2 )
30    continue
c
      return
      end
