      parameter( nr =2000    , ! nr must be even
     x           nz =1001    ,
     x           nx =201     ,
     x           ng =51      ,
     x           ne =15      ,
     x           nlm=18      ,
     x           ncf=2       ,
     x           hr =0.01d0  ,
     x           xx0=1.25d0  ,
c
     x           xf=200.d0         ,
     x           pi =dacos(-1.d0 ) , 
     x           pi2=pi*pi         ,
     x           pi3=pi*pi2        ,
c
     x           hbarc=197.32697d0                ,
     x           ampi0=134.9766d0/hbarc           ,
     x           ampic=139.5702d0/hbarc           ,
     x           amp  =( ampi0+2.d0*ampic )/3.d0  ,
     x           amd  =293.1d0/hbarc              ,
c
     x           ap02=ampi0**2 ,
     x           apc2=ampic**2 ,
     x           amp2=amp**2   ,
     x           amd2=amd**2   ,
c
     x           gan0=  1.29d0       ,
     x           gan2=gan0**2        ,
     x           gan4=gan2**2        ,
     x           fpp0=184.80d0/hbarc ,
c
     x           ct2 =( gan0/fpp0 )**2     ,
     x           ct4 =ct2**2/( 24.d0*pi2 ) ,
     x           ct48=ct2**2/(  8.d0*pi3)  ,
c
     x           dm2=- 5.d0+ 4.d0/gan0**2+1.d0/gan0**4 ,
     x           dk2=-23.d0+10.d0/gan0**2+1.d0/gan0**4 ,
c  
     x           aftt =ct2/(12.d0*pi) ,
     x           aft  =ct48*amp       ,
c
     x           afcst=aftt*amp2  ,
     x           afcs =4.d0*aft   ,
     x           afct=aft/gan0**4 )
