c    <Copyright (c) 1999 by Mahendra Bhagwat>
c
      subroutine avgiflow(nr,np,nw,pw,nz,ns,nk,dp,dz,ip,anb,crd,rad,
     &                 dcy,kinv,om,lin,gav,gv,rc0,asr,b0r,b1c,b1s,
     &                 rr_x,rr_y,rr_z,pcx,pcy,pcz,cpx,cpy,cpz,
     &                 sco,sce,gsc,csc,compress,lengthnw,
     &                 pif,zif,bvx,bvy,bvz,rcb,gb,BICM,NWICM,
     &                 eqtim)
c
c     Time averaged induced velocity field in the TPP
c----------------------------------------------------------------------
c
      IMPLICIT none
c
        include 'adimpar.inc'
        include 'kdimpar.inc'
        include 'rdimpar.inc'
        include 'sdimpar.inc'
        include 'wdimpar.inc'
        include 'zdimpar.inc'
c
      INTEGER nr,np,nw,pw,nz,ns,nk,r,p,s,pn,dnp,
     &        r0,p0,s0,r1,p1,w1,z1,npl
c
      DOUBLE PRECISION dp,dz,ip(rdim),anb,dcy,kinv,om,lin,crd,rad,
     &       gav,trb,az,sn,vbx,vby,vbz,pi,
     &       t11,t12,t13,
     &       t21,t22,t23,
     &       t31,t32,t33,
     &       sce,b0
c
      DOUBLE PRECISION asr(rdim),b0r(rdim),rr_x(rdim),rr_y(rdim),
     &                                     rr_z(rdim)
c
      DOUBLE PRECISION b1c(rdim),b1s(rdim)
c
      DOUBLE PRECISION cpx(rdim,adim,sdim),
     &       cpy(rdim,adim,sdim),
     &       cpz(rdim,adim,sdim),
     &       vbtx(rdim,adim,sdim),
     &       vbty(rdim,adim,sdim),
     &       vbtz(rdim,adim,sdim),
     &       vbtx1(rdim,adim,sdim),
     &       vbty1(rdim,adim,sdim),
     &       vbtz1(rdim,adim,sdim),
     &       vibx(rdim,adim,sdim),
     &       viby(rdim,adim,sdim),
     &       vibz(rdim,adim,sdim),
     &       vibx1(rdim,adim,sdim),
     &       viby1(rdim,adim,sdim),
     &       vibz1(rdim,adim,sdim)

c
c   instantaneour velocities over the disk
c   along the Psi=0 line
      DOUBLE PRECISION vbtxi(rdim,adim,sdim),
     &       vbtyi(rdim,adim,sdim),
     &       vbtzi(rdim,adim,sdim),
     &       vibxi(rdim,adim,sdim),
     &       vibyi(rdim,adim,sdim),
     &       vibzi(rdim,adim,sdim)
c
c
      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim),
     &       rc0(rdim,adim,wdim,kdim),
     &       pcx(rdim,adim,wdim,zdim),
     &       pcy(rdim,adim,wdim,zdim),
     &       pcz(rdim,adim,wdim,zdim),
     &       eqtim(rdim,adim,wdim,zdim,kdim)
c
      CHARACTER*1 sco
c
      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim),
     &       csc(rdim,adim,wdim,zdim,kdim)
c
      INTEGER lr,lp,lk
c
      INTEGER pif,zif,Nb
c
      DOUBLE PRECISION bvx(rdim,adim,sdim),
     &       bvy(rdim,adim,sdim),
     &       bvz(rdim,adim,sdim)
      DOUBLE PRECISION rcb
      DOUBLE PRECISION gb(rdim,adim,sdim)
      INTEGER s1,firstz,dp90,p270
      DOUBLE PRECISION NWICM(sdim,sdim), BICM(sdim,sdim)
c
      INTEGER peff,blade
c
      INTEGER lengthnw(sdim)
c
      DOUBLE PRECISION planf
c
      DOUBLE PRECISION compress(rdim,adim,sdim)
c
      DOUBLE PRECISION b2c(rdim),b2s(rdim),b3c(rdim),b3s(rdim),
     &                 b4c(rdim),b4s(rdim)
c
      DOUBLE PRECISION factor,rc,rc_2p,rc_4p,zet0,wage
      INTEGER k
c
      external planf
c
c----------------------------------------------------------------------
c
      common/pidata/pi
      save/pidata/
c
      common/hhflap/b2c,b2s,b3c,b3s,b4c,b4s
      save/hhflap/
c----------------------------------------------------------------------
c
      dp90=nint(pi/(2.D+0*(dp/DBLE(pif))))
      p270=3*dp90+1
c
c     Initialize total induced velocity at every blade
c     control point to zero so it isn't re-summed:
      do r0=1,nr
      do p0=1,np
      do s0=1,ns

            vbtx(r0,p0,s0)=0.D+0
            vbty(r0,p0,s0)=0.D+0
            vbtz(r0,p0,s0)=0.D+0
            vbtx1(r0,p0,s0)=0.D+0
            vbty1(r0,p0,s0)=0.D+0
            vbtz1(r0,p0,s0)=0.D+0
            vbtxi(r0,p0,s0)=0.D+0
            vbtyi(r0,p0,s0)=0.D+0
            vbtzi(r0,p0,s0)=0.D+0

      end do
      end do
      end do
c----------------------------------------------------------------------
c
      dnp=nint(anb*100000.D+0/(dp/DBLE(pif))/100000.D+0)
c
      Nb=nint(2*pi*100000.D+0/anb/100000.D+0)
c
c----------------------------------------------------------------------
c
c     For peace of mind, initialize:
      vbx=0.D+0
      vby=0.D+0
      vbz=0.D+0
c----------------------------------------------------------------------
c     Vortex/blade interactions:
          gav=0.D+0
          do lr=1,nr
          do lp=1,np
             do lk=1,nk
             gav=gav+abs(gv(lr,lp,1,1,lk))
             end do
          end do
          end do
          gav=gav/DBLE(nr*np)

c-mb re-compute Gamma_avg and therefore Delta
          trb=1.D+0+1.D+0*dcy*gav/kinv

c
c     Influenced indices:
      do r0=1,nr
      do p0=1,np   !p270,dp90 
      npl=p0+np-dnp
      do s0=1,ns

c        Forcing indices:
         do r1=1,nr

         sn=(-1.D+0)**(r1-1)

         do pn=1,np,1
c-mb note this is time-averaged 
         p1=pn
         firstz=2
c-mb         do w1=1,nw    !full-span wake
         do w1=1,nw-pw !no inboard wake
         if(p1.eq.p0) firstz=2+lengthnw(nw+1-w1)
c-mb this is because at the blade TV is accounted for in the near-wake
         do z1=firstz,nz,1

          factor=1.D+0
          k=1
c>> mb more cores?
          wage=DBLE(z1-1)*dz/DBLE(zif)
          trb=1.D+0+1.D+0*dcy*factor*abs(gv(r1,p1,w1,z1,1)/kinv)
          trb=1.D+0+1.D+0*dcy*factor*gav/kinv
          zet0=(rc0(r1,p1,w1,k)/.00855D+0)**2*om/trb/Nb
          rc=0.00855D+0*sqrt(trb*(wage+zet0)*Nb/om)

c-S Modified core radius approach...
          rc=0.00855D+0*sqrt(trb*eqtim(r1,p1,w1,z1,1)*Nb/om)


            call vbind(r0,p0,s0,r1,p1,w1,z1,np,nz,nk,dz,
     &                 gv,rc,lin,trb,om,rc0,cpx,cpy,cpz,
     &                 pcx,pcy,pcz,vbx,vby,vbz,zif,Nb)

c           Sum velocity contributions for total
c           induced velocity at each control point:
            vbtx(r0,p0,s0)=vbtx(r0,p0,s0)+vbx*sn
            vbty(r0,p0,s0)=vbty(r0,p0,s0)+vby*sn
            vbtz(r0,p0,s0)=vbtz(r0,p0,s0)+vbz*sn
            vbtx1(r0,p0,s0)=vbtx1(r0,p0,s0)+vbx*sn
            vbty1(r0,p0,s0)=vbty1(r0,p0,s0)+vby*sn
            vbtz1(r0,p0,s0)=vbtz1(r0,p0,s0)+vbz*sn
c-mb instantaneous velocities
            if(p0.eq.1) then
            do blade=1,Nb,1
             peff=p1+(blade-1)*dnp
             if(peff.gt.np) peff=peff-np
             vbtxi(r0,peff,s0)=vbtxi(r0,peff,s0)+vbx*sn
             vbtyi(r0,peff,s0)=vbtyi(r0,peff,s0)+vby*sn
             vbtzi(r0,peff,s0)=vbtzi(r0,peff,s0)+vbz*sn
            end do
            endif

         end do !z1
         end do !w1
         end do !p1
         end do !r1

      end do
      end do
      end do
c----------------------------------------------------------------------
c
c     For peace of mind, reinitialize:
      vbx=0.D+0
      vby=0.D+0
      vbz=0.D+0
c----------------------------------------------------------------------
c
c     Shed circulation effects:

      if (sco .eq. 'y') then
      if (nw  .gt.  1 ) then

c     Influenced indices:
      do r0=1,nr
      do p0=1,np  !p270,dp90  
      npl=p0+np-dnp
      do s0=1,ns

c        Forcing indices:
         do r1=1,nr
         do pn=1,np,1
         p1=pn
         firstz=2
c-mb         do w1=2,nw    !full-span wake
         do w1=2,nw-pw !no inboard wake
         if(p1.eq.p0) firstz=2+lengthnw(nw+1-w1)
c-mb this is because at the blade TV is accounted for in the near-wake
         do z1=firstz,nz,1


          factor=1.D+0
          k=1
c>> mb more cores?
          wage=DBLE(z1-1)*dz/DBLE(zif)
          trb=1.D+0+1.D+0*dcy*factor*abs(gv(r1,p1,w1,z1,1)/kinv)
          trb=1.D+0+1.D+0*dcy*factor*gav/kinv
          zet0=(rc0(r1,p1,w1,k)/.00855D+0)**2*om/trb/Nb
          rc=0.00855D+0*sqrt(trb*(wage+zet0)*Nb/om)

c-S Modified core radius approach...
          rc=0.00855D+0*sqrt(trb*eqtim(r1,p1,w1,z1,1)*Nb/om)

            call sbind(r0,p0,s0,r1,p1,w1,z1,np,nk,dz,om,rc0,
     &                 gv,rc,lin,trb,sce,gsc,csc,cpx,cpy,cpz,
     &                 pcx,pcy,pcz,vbx,vby,vbz,zif)

c           Sum velocity contributions for total induced
c           velocity at each collocation point:
            vbtx(r0,p0,s0)=vbtx(r0,p0,s0)+vbx
            vbty(r0,p0,s0)=vbty(r0,p0,s0)+vby
            vbtz(r0,p0,s0)=vbtz(r0,p0,s0)+vbz
            vbtx1(r0,p0,s0)=vbtx1(r0,p0,s0)+vbx
            vbty1(r0,p0,s0)=vbty1(r0,p0,s0)+vby
            vbtz1(r0,p0,s0)=vbtz1(r0,p0,s0)+vbz
c-mb instantaneous velocities
            if(p0.eq.1) then
            do blade=1,Nb,1
             peff=p1+(blade-1)*dnp
             if(peff.gt.np) peff=peff-np
             vbtxi(r0,peff,s0)=vbtxi(r0,peff,s0)+vbx
             vbtyi(r0,peff,s0)=vbtyi(r0,peff,s0)+vby
             vbtzi(r0,peff,s0)=vbtzi(r0,peff,s0)+vbz
            end do
            endif

         end do
         end do
         end do
         end do

      end do
      end do
      end do

      endif
      endif
c----------------------------------------------------------------------
c
c     For Dr. Bugi's peace of mind, reinitialize:
      vbx=0.D+0
      vby=0.D+0
      vbz=0.D+0
c----------------------------------------------------------------------
c
c     Induced inflow : velocity induced by the near-wake
c
c     Influenced indices:
      do r0=1,nr
      do p0=1,np   !p270,dp90  
      npl=p0+np-dnp
      do s0=1,ns

c        Forcing indices:
         do r1=r0,r0
         do pn=p0,p0
c-mb note that this contribution is only from one blade
c-mb lateron I include Nb
         p1=pn
         if (p1 .gt. np) p1=p1-np
         do s1=1,ns

c-mb near-wake is assumed to induce only downwash 
c-mb     (no chord- / span-wise)
            vbx=0.D+0
            vby=0.D+0
            vbz=-NWICM(s0,s1)*gb(r1,p1,s1)
c-mb note the negative sign

c-mb also, this is relative to TPP ( tilt by shaft angle :)
            vbx=vbz*sin(-asr(r1))
            vbz=vbz*cos(-asr(r1))

c-mb instantaneous velocities
            if(p0.eq.1) then
            do blade=1,Nb,1
             peff=p1+(blade-1)*dnp
             if(peff.gt.np) peff=peff-np
             vbtxi(r0,peff,s0)=vbtxi(r0,peff,s0)+vbx
             vbtyi(r0,peff,s0)=vbtyi(r0,peff,s0)+vby
             vbtzi(r0,peff,s0)=vbtzi(r0,peff,s0)+vbz
            end do
            endif

c-mb  assume this to be over the blade area! 
c-mb     Chord*DeltaS/(2*pi*cpx*DeltaS)
c-mb            vbz=vbz*DBLE(Nb)*crd*planf(cpx(r0,1,s0)/rad,taperst,taper)/
c-mb     &          (2.D+0*pi*cpx(r0,1,s0))
c-mb  this is not necessary as in effect what I am doing is replacing
c-mb  first free segment with the near-wake :)

c           Sum velocity contributions for total induced
c           velocity at each collocation point:
            vbtx(r0,p0,s0)=vbtx(r0,p0,s0)+vbx
            vbty(r0,p0,s0)=vbty(r0,p0,s0)+vby
            vbtz(r0,p0,s0)=vbtz(r0,p0,s0)+vbz
            vbtx1(r0,p0,s0)=vbtx1(r0,p0,s0)+vbx
            vbty1(r0,p0,s0)=vbty1(r0,p0,s0)+vby
            vbtz1(r0,p0,s0)=vbtz1(r0,p0,s0)+vbz

         end do
         end do
         end do

      end do
      end do
      end do

c----------------------------------------------------------------------
c-mb till here the velocity was summed over the azimuth so divide 
c-mb  and multiply by number of blades

      do r0=1,nr
      do s0=1,ns
      do p0=1,np    !p270,dp90  
            vbtx(r0,p0,s0)=vbtx(r0,p0,s0)*DBLE(Nb)/DBLE(np)
            vbty(r0,p0,s0)=vbty(r0,p0,s0)*DBLE(Nb)/DBLE(np)
            vbtz(r0,p0,s0)=vbtz(r0,p0,s0)*DBLE(Nb)/DBLE(np)
            vbtx1(r0,p0,s0)=vbtx1(r0,p0,s0)*DBLE(Nb)/DBLE(np)
            vbty1(r0,p0,s0)=vbty1(r0,p0,s0)*DBLE(Nb)/DBLE(np)
            vbtz1(r0,p0,s0)=vbtz1(r0,p0,s0)*DBLE(Nb)/DBLE(np)
      end do
      end do
      end do


c----------------------------------------------------------------------
c
c     For peace of mind, initialize:
      vbx=0.D+0
      vby=0.D+0
      vbz=0.D+0
c----------------------------------------------------------------------
c
c     Bound-vortex induced inflow (momentum addition)
c-mb velocity induced by the bound vortices: this is 
c-mb  because this is a momentum change
c-mb  we want the velocity induced by the bound vortices at the tip
c-mb  and at the 3/4 chord point (where we have the evaluation pts)
c-mb  note that we want the velocity induced by all OTHER segments!
c
c     Influenced indices:
      do r0=1,nr
      do p0=1,np   !p270,dp90  
      npl=p0+np-dnp
      do s0=1,ns

c     Forcing indices:
         do r1=r0,r0
         do p1=p0,p0,1
         do s1=1,ns

c-mb velocity induced by the bound vortex segments
            vbx=0.D+0
            vby=0.D+0
            vbz=-BICM(s0,s1)*gb(r1,p1,s1)
c-mb note the negative sign

c-mb also, this is relative to TPP ( tilt by shaft angle :)
            vbx=vbz*sin(-asr(r1))
            vbz=vbz*cos(-asr(r1))

c-mb instantaneous velocities
            if(p0.eq.1) then
            do blade=1,Nb,1
             peff=p1+(blade-1)*dnp
             if(peff.gt.np) peff=peff-np
             vbtxi(r0,peff,s0)=vbtxi(r0,peff,s0)+vbx
             vbtyi(r0,peff,s0)=vbtyi(r0,peff,s0)+vby
             vbtzi(r0,peff,s0)=vbtzi(r0,peff,s0)+vbz
            end do
            endif

c-mb  assume this to be over the blade area! 
c-mb        Chord*DeltaS/(2*pi*cpx*DeltaS)
c-mb   
            vbz=vbz*DBLE(Nb)/(2.D+0*pi)

c-mb ok, now I think I dont want Nb
            vbz=vbz/DBLE(Nb)  

c           Sum velocity contributions for total induced
c           velocity at each collocation point:
            vbtx(r0,p0,s0)=vbtx(r0,p0,s0)+vbx
            vbty(r0,p0,s0)=vbty(r0,p0,s0)+vby
            vbtz(r0,p0,s0)=vbtz(r0,p0,s0)+vbz


         end do
         end do
         end do


      end do
      end do
      end do



c----------------------------------------------------------------------
c
c     Fixed to rotating frame transform:
      do r=1,nr
      sn=(-1.D+0)**(r-1)
      do p=1,np

         az=(ip(r)+DBLE(p-1)*dp) - ip(r)
c-mb VFF correctn
         az=(ip(r)+DBLE(p-1)*dp/DBLE(pif)) - ip(r)
         b0=b0r(r)  ! bagai had this for whatever reasons in rindv.f
         b0=b0r(r)+b1c(r)*cos(az)+b1s(r)*sin(az)
     &            +b2c(r)*cos(2.D+0*az)+b2s(r)*sin(2.D+0*az)
     &            +b3c(r)*cos(3.D+0*az)+b3s(r)*sin(3.D+0*az)
     &            +b4c(r)*cos(4.D+0*az)+b4s(r)*sin(4.D+0*az)

         az=sn*(az + ip(r))
c-mb sign of rotation for coord transform

         t11= cos(-asr(r))*cos(az)*cos(b0)+
     &        sin(-asr(r))*sin(b0)
         t12= sin(az)*cos(b0)
         t13=-sin(-asr(r))*cos(az)*cos(b0)+
     &        cos(-asr(r))*sin(b0)

         t21=-cos(-asr(r))*sin(az)*sn
         t22= cos(az)*sn
         t23= sin(-asr(r))*sin(az)*sn

         t31=-cos(-asr(r))*cos(az)*sin(b0)+
     &        sin(-asr(r))*cos(b0)
         t32=-sin(az)*sin(b0)
         t33= sin(-asr(r))*cos(az)*sin(b0)+
     &        cos(-asr(r))*cos(b0)

c----------------------------------------------------------------------|
      do s=1,ns
       vibx(r,p,s)=   t11*vbtx(r,p,s)+t12*vbty(r,p,s)+t13*vbtz(r,p,s)
       viby(r,p,s)=+( t21*vbtx(r,p,s)+t22*vbty(r,p,s)+t23*vbtz(r,p,s) )
       vibz(r,p,s)=   t31*vbtx(r,p,s)+t32*vbty(r,p,s)+t33*vbtz(r,p,s)

       vibx1(r,p,s)=
     &   t11*vbtx1(r,p,s)+t12*vbty1(r,p,s)+t13*vbtz1(r,p,s)
       viby1(r,p,s)=
     &+( t21*vbtx1(r,p,s)+t22*vbty1(r,p,s)+t23*vbtz1(r,p,s) )
       vibz1(r,p,s)=
     &   t31*vbtx1(r,p,s)+t32*vbty1(r,p,s)+t33*vbtz1(r,p,s)

       vibxi(r,p,s)=
     &   t11*vbtxi(r,p,s)+t12*vbtyi(r,p,s)+t13*vbtzi(r,p,s)
       vibyi(r,p,s)=
     &+( t21*vbtxi(r,p,s)+t22*vbtyi(r,p,s)+t23*vbtzi(r,p,s) )
       vibzi(r,p,s)=
     &   t31*vbtxi(r,p,s)+t32*vbtyi(r,p,s)+t33*vbtzi(r,p,s)
      end do
c----------------------------------------------------------------------|

      end do
      end do

      open(71,file='VBXAVG.dat',status='unknown')
      open(72,file='VBYAVG.dat',status='unknown')
      open(73,file='VBZAVG.dat',status='unknown')

          do r=1,nr
 
          do p=1,1
          do s=ns,1,-1
             write(71,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
     &                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)),
     &                  vibx(r,p,s)
             write(72,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
     &                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)),
     &                  viby(r,p,s)
             write(73,'(3(E15.7,2x))')
     &                  (cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
     &                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)),
     &                  vibz(r,p,s),vibz1(r,p,s)
          end do
          end do

             write(71,*)
             write(72,*)
             write(73,*)

          do p=1+int(pi/(dp/DBLE(pif))),1+int(pi/(dp/DBLE(pif)))
c-mb VFF correctn          ^^^^^^^^^^               ^^^^^^^^^^
          do s=1,ns
             write(71,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
     &                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)),
     &                  vibx(r,p,s)
             write(72,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
     &                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)),
     &                  viby(r,p,s)
             write(73,'(3(E15.7,2x))')
     &                  (cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
     &                 -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)),
     &                  vibz(r,p,s),vibz1(r,p,s)
          end do
          end do

             write(71,*)
             write(72,*)
             write(73,*)

c-mb VFF correctn              __________
          do p=int(pi/(2.D+0*(dp/DBLE(pif))))+1,
     &         int(pi/(2.D+0*(dp/DBLE(pif))))+1
c-mb VFF correctn              ^^^^^^^^^^
c----------------------------------------------------------------------|
          do s=ns,1,-1
             write(71,*)(cpy(r,p,s)-rr_y(r))/rad,
     &                  vibx(r,p,s)
             write(72,*)(cpy(r,p,s)-rr_y(r))/rad,
     &                  viby(r,p,s)
           write(73,'(3(E15.7,2x))')(cpy(r,p,s)-rr_y(r))/rad,
     &                  vibz(r,p,s),vibz1(r,p,s)
          end do
          end do

             write(71,*)
             write(72,*)
             write(73,*)

c-mb VFF correctn              __________                __________
          do p=int(pi/(2.D+0*(dp/DBLE(pif))))+1+int(pi/(dp/DBLE(pif))),
     &         int(pi/(2.D+0*(dp/DBLE(pif))))+1+int(pi/(dp/DBLE(pif)))
c-mb VFF correctn              ^^^^^^^^^^                ^^^^^^^^^^
          do s=1,ns
             write(71,*)(cpy(r,p,s)-rr_y(r))/rad,
     &                  vibx(r,p,s)
             write(72,*)(cpy(r,p,s)-rr_y(r))/rad,
     &                  viby(r,p,s)
           write(73,'(3(E15.7,2x))')(cpy(r,p,s)-rr_y(r))/rad,
     &                  vibz(r,p,s),vibz1(r,p,s)
          end do
          end do

          end do

      close(71)
      close(72)
      close(73)
c----------------------------------------------------------------------|
c-mb  also write out the avg inflow w/o momentum addition may be needed
      open(74,file='VBZAVG_PSI.dat',status='unknown')
      do r=1,nr
      do p=1,np
      do s=1,ns
      write(74,174)cpx(r,p,s),cpy(r,p,s),cpz(r,p,s)
     &            ,vibx(r,p,s),viby(r,p,s),vibz(r,p,s),vibz1(r,p,s)
      end do
      write(74,*)
      write(74,*)
      end do
      write(74,*)
      end do
      close(74)
174   format(7(E12.6,1X))
c----------------------------------------------------------------------|
      open(75,file='VBINST.dat',status='unknown')
      do r=1,nr
      do p=1,np
      do s=1,ns
      write(75,175)cpx(r,p,s),cpy(r,p,s),cpz(r,p,s)
     &            ,vibxi(r,p,s),vibyi(r,p,s),vibzi(r,p,s)
      end do
      write(75,*)
      end do
      write(75,*)
      end do
      close(75)
175   format(6(E13.6,2X))
c----------------------------------------------------------------------|
c
      write(*,*)'finished time-averaged velos'
      return
      end 
