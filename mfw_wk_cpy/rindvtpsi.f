C-$	$Id: rindvtpsi.f,v 1.5.2.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine rindvtpsi(nr,np,nw,pw,nz,ns,nk,dp,dz,ip,anb,crd,rad,
     $     dcy,kinv,om,lin,gav,gv,rc0,asr,beta, rr_x,rr_y,rr_z,pcx,pcy
     $     ,pcz,cpx,cpy,cpz, sco,sce,gsc,csc,vibx,viby,vibz, vibxi,vibyi
     $     ,vibzi,compress,lengthnw, pif,zif,bvx,bvy,bvz,rcb,gb,BICM
     $     ,NWICM,p0, eqtim,lsr)


c     Induced velocity field at blade contol points.

      IMPLICIT none

      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER nr,np,nw,pw,nz,ns,nk,r,p,s,pn,dnp, r0,p0,s0,r1,p1,w1,z1
     $     ,npl

      DOUBLE PRECISION dp,dz,ip(rdim),anb,dcy,kinv,om,lin,crd,rad, gav
     $     ,trb,az,sn,vbx,vby,vbz,pi, t11,t12,t13, t21,t22,t23, t31,t32
     $     ,t33, sce,b0

      DOUBLE PRECISION asr(rdim),beta(rdim,adim),rr_x(rdim),rr_y(rdim),
     $     rr_z(rdim),lsr(rdim)

      DOUBLE PRECISION cpx(rdim,adim,sdim), cpy(rdim,adim,sdim),
     $     cpz(rdim,adim,sdim), vbtx(rdim,adim,sdim), vbty(rdim,adim
     $     ,sdim), vbtz(rdim,adim,sdim), vbtxin(rdim,adim,sdim),
     $     vbtyin(rdim,adim,sdim), vbtzin(rdim,adim,sdim), vibx(rdim
     $     ,adim,sdim), viby(rdim,adim,sdim), vibz(rdim,adim,sdim),
     $     vibxi(rdim,adim,sdim), vibyi(rdim,adim,sdim), vibzi(rdim,adim
     $     ,sdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), rc0(rdim,adim,wdim
     $     ,kdim), pcx(rdim,adim,wdim,zdim), pcy(rdim,adim,wdim,zdim),
     $     pcz(rdim,adim,wdim,zdim), eqtim(rdim,adim,wdim,zdim,kdim)

      CHARACTER*1 sco

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim), csc(rdim,adim,wdim
     $     ,zdim,kdim)

      INTEGER lr,lp,lk

      INTEGER pif,zif,Nb

      DOUBLE PRECISION bvx(rdim,adim,sdim), bvy(rdim,adim,sdim),
     $     bvz(rdim,adim,sdim)
      DOUBLE PRECISION rcb
      DOUBLE PRECISION gb(rdim,adim,sdim)
      INTEGER s1
      DOUBLE PRECISION NWICM(sdim,sdim), BICM(sdim,sdim)

      INTEGER zstart,zend,zstep

      DOUBLE PRECISION compress(rdim,adim,sdim)

      DOUBLE PRECISION trbin

      DOUBLE PRECISION sumb,sumnw
      INTEGER i,j

      INTEGER lengthnw(sdim)

      DOUBLE PRECISION factor,rc,rc_2p,rc_4p,zet0,wage
      INTEGER k
      INTEGER rotgeo


      common/pidata/pi
      save/pidata/

      common/rotdata/rotgeo
      save/rotdata/


      dnp=nint(anb/(dp/DBLE(pif)))

      Nb=nint(2*pi/anb)

c     Initialize total induced velocity at every blade
c     control point to zero so it isn't re-summed:
      do r0=1,nr
         do s0=1,ns

            vbtx(r0,p0,s0)=0.D+0
            vbty(r0,p0,s0)=0.D+0
            vbtz(r0,p0,s0)=0.D+0
            vbtxin(r0,p0,s0)=0.D+0
            vbtyin(r0,p0,s0)=0.D+0
            vbtzin(r0,p0,s0)=0.D+0

         end do
      end do

c     For peace of mind, initialize:
      vbx=0.D+0
      vby=0.D+0
      vbz=0.D+0

c     Vortex/blade interactions:
c-mb  only for the tip vortex
      gav=0.D+0
      do lr=1,nr
         do lp=1,np
            do lk=1,nk
               gav=gav+abs(gv(lr,lp,1,1,lk))
            end do
         end do
      end do
      gav=gav/DBLE(nr*np)

      trb=1.D+0+dcy*gav/kinv

c     Influenced indices:
      do r0=1,nr
         npl=p0+np-dnp
         do s0=1,ns

c     Forcing indices:
            do r1=1,nr

               sn=1.D+0
               if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r1-1)

               do pn=p0,npl,dnp
                  p1=pn
                  if (p1 .gt. np) p1=p1-np
                  do w1=1,nw-pw !no inboard wake
                     zstart=2
                     zend=nz
                     zstep=1
                     if(nz .gt. (nint(2.D+0*2.D+0*pi/dz)*zif+1)) zend
     $                    =(nint(2.D+0*2.D+0*pi/dz)*zif+1)
                     zstart=3!+lengthnw(nw+1-w1)
                     zend=nz

                     do z1=zstart,zend,zstep
                        k=1
                        rc=0.00855D+0*sqrt(trb*eqtim(r1,p1,w1,z1,k)*Nb
     $                       /om)

                        call vbind(r0,p0,s0,r1,p1,w1,z1,np,nz,nk,dz, gv
     $                       ,rc,lin,trbin,om,rc0,cpx,cpy,cpz, pcx,pcy
     $                       ,pcz,vbx,vby,vbz,zif,Nb)

c     Sum velocity contributions for total
c     induced velocity at each control point:
                        vbtx(r0,p0,s0)=vbtx(r0,p0,s0)+vbx*sn
                        vbty(r0,p0,s0)=vbty(r0,p0,s0)+vby*sn
                        vbtz(r0,p0,s0)=vbtz(r0,p0,s0)+vbz*sn
                        vbtxin(r0,p0,s0)=vbtxin(r0,p0,s0)+vbx*sn
                        vbtyin(r0,p0,s0)=vbtyin(r0,p0,s0)+vby*sn
                        vbtzin(r0,p0,s0)=vbtzin(r0,p0,s0)+vbz*sn

                     end do
                  end do
               end do
            end do
         end do
      end do

c     For peace of mind, reinitialize:
      vbx=0.D+0
      vby=0.D+0
      vbz=0.D+0
c----------------------------------------------------------------------

c     Shed circulation effects:
      if (sco .eq. 'y') then
         if (nw-pw  .gt.  1 ) then !no inboard wake

c     Influenced indices:
            do r0=1,nr
               npl=p0+np-dnp
               do s0=1,ns

c     Forcing indices:
                  do r1=1,nr
                     do pn=p0,npl,dnp 
                        p1=pn
                        if (p1 .gt. np) p1=p1-np
                        do w1=2,nw-pw !no inboard wake
                           zstart=3
                           zend=nz
                           zstep=1
                           zstart=2+lengthnw(nw+1-w1)

                           if(nz .gt. (nint(2.D+0*2.D+0*pi/dz)*zif+1))
     $                          zend=(nint(2.D+0*2.D+0*pi/dz)*zif+1)

                           do z1=zstart,zend,zstep
                              k=1
                              rc=0.00855D+0*sqrt(trb*eqtim(r1,p1,w1,z1,k
     $                             )*Nb/om)

                              call sbind(r0,p0,s0,r1,p1,w1,z1,np,nk,dz
     $                             ,om,rc0,gv,rc,lin,trbin,sce,gsc,csc
     $                             ,cpx,cpy,cpz,pcx,pcy,pcz,vbx,vby,vbz
     $                             ,zif)

c     Sum velocity contributions for total induced
c     velocity at each collocation point:
                              vbtx(r0,p0,s0)=vbtx(r0,p0,s0)+vbx
                              vbty(r0,p0,s0)=vbty(r0,p0,s0)+vby
                              vbtz(r0,p0,s0)=vbtz(r0,p0,s0)+vbz
                              vbtxin(r0,p0,s0)=vbtxin(r0,p0,s0)+vbx
                              vbtyin(r0,p0,s0)=vbtyin(r0,p0,s0)+vby
                              vbtzin(r0,p0,s0)=vbtzin(r0,p0,s0)+vbz

                           end do
                        end do
                     end do
                  end do

               end do
            end do

         endif
      endif

      vbx=0.D+0
      vby=0.D+0
      vbz=0.D+0

c     Induced inflow : velocity induced by the near-wake
c      this needs to be included for calculating AOAinduced
c      otherwise, W-L solver has this inherently coupled.

c     Influenced indices:
      do r0=1,nr
         npl=p0+np-dnp
         do s0=1,ns

c     Forcing indices:
            do r1=r0,r0
               do pn=p0,p0
                  p1=pn
                  if (p1 .gt. np) p1=p1-np
                  do s1=1,ns

                     vbx=0.D+0
                     vby=0.D+0
                     vbz=-(NWICM(s0,s1))*gb(r1,p1,s1)
     $                    *compress(r1,p1,s1)
                     vbx=vbz*sin(-asr(r1))*cos(lsr(r1))
                     vbz=vbz*cos(-asr(r1))*cos(lsr(r1))

c     Sum velocity contributions for total induced
c     velocity at each collocation point:
                     vbtxin(r0,p0,s0)=vbtxin(r0,p0,s0)+vbx
                     vbtyin(r0,p0,s0)=vbtyin(r0,p0,s0)+vby
                     vbtzin(r0,p0,s0)=vbtzin(r0,p0,s0)+vbz
                  end do        !s1
               end do           !pn/p1
            end do              !r1

         end do                 !s0
      end do                    !r0


c     Fixed to rotating frame transform:
      do r=1,nr
         p=p0
         sn=1.D+0
         if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r-1)

         az=(DBLE(p-1)*dp/DBLE(pif))
         b0=beta(r,p0)

         az=sn*(az + ip(r))

C        t11= cos(-asr(r))*cos(az)*cos(b0)+ sin(-asr(r))*sin(b0)
C        t12= sin(az)*cos(b0)
C        t13=-sin(-asr(r))*cos(az)*cos(b0)+ cos(-asr(r))*sin(b0)
C
C        t21=-cos(-asr(r))*sin(az)*sn
C        t22= cos(az)*sn
C        t23= sin(-asr(r))*sin(az)*sn
C
C        t31=-cos(-asr(r))*cos(az)*sin(b0)+ sin(-asr(r))*cos(b0)
C        t32=-sin(az)*sin(b0)
C        t33= sin(-asr(r))*cos(az)*sin(b0)+ cos(-asr(r))*cos(b0)

C/MR  Adding lateral shaft tilt angle          
         t11= cos(-asr(r))*cos(az)*cos(b0)
     $          -sn*sin(-asr(r))*sin(lsr(r))*sin(az)*cos(b0) 
     $          +sin(-asr(r))*cos(lsr(r))*sin(b0)
         t12= sn*cos(lsr(r))*sin(az)*cos(b0)+sin(lsr(r))*sin(b0)
         t13=-sin(-asr(r))*cos(az)*cos(b0)
     $          -sn*cos(-asr(r))*sin(lsr(r))*sin(az)*cos(b0)
     $          +cos(-asr(r))*cos(lsr(r))*sin(b0)

         t21=-cos(-asr(r))*sin(az)-sn*sin(asr(r))*sin(lsr(r))*cos(az)
         t22= cos(lsr(r))*cos(az)*sn
         t23= sin(-asr(r))*sin(az)-sn*cos(asr(r))*sin(lsr(r))*cos(az)

         t31=-cos(-asr(r))*cos(az)*sin(b0)
     $          +sn*sin(-asr(r))*sin(lsr(r))*sin(az)*sin(b0) 
     $          +sin(-asr(r))*cos(lsr(r))*cos(b0)
         t32=-sn*cos(lsr(r))*sin(az)*sin(b0)+sin(lsr(r))*cos(b0)
         t33= sin(-asr(r))*cos(az)*sin(b0)
     $          +sn*cos(asr(r))*sin(lsr(r))*sin(az)*sin(b0) 
     $          +cos(-asr(r))*cos(lsr(r))*cos(b0)

         do s=1,ns

            vibx(r,p,s)=   t11*vbtx(r,p,s)+t12*vbty(r,p,s)+t13*vbtz(r,p
     $           ,s)
            viby(r,p,s)=+( t21*vbtx(r,p,s)+t22*vbty(r,p,s)+t23*vbtz(r,p
     $           ,s) )
            vibz(r,p,s)=   t31*vbtx(r,p,s)+t32*vbty(r,p,s)+t33*vbtz(r,p
     $           ,s)
            vibxi(r,p,s)=   t11*vbtxin(r,p,s) +t12*vbtyin(r,p,s) +t13
     $           *vbtzin(r,p,s)
            vibyi(r,p,s)=+( t21*vbtxin(r,p,s) +t22*vbtyin(r,p,s) +t23
     $           *vbtzin(r,p,s) )
            vibzi(r,p,s)=   t31*vbtxin(r,p,s) +t32*vbtyin(r,p,s) +t33
     $           *vbtzin(r,p,s)
         end do                 !s

      end do                    !r

      return
      end 
