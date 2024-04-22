C-$	$Id: rindvt.f,v 1.5.2.2 2005/06/11 15:57:37 shreyas Exp $	

      subroutine rindvt(nr,np,nw,pw,nz,ns,nk,dp,dz,ip,anb,crd,rad, dcy
     $     ,kinv,om,lin,gav,gv,rc0,asr,b0r,b1c,b1s, rr_x,rr_y,rr_z,pcx
     $     ,pcy,pcz,cpx,cpy,cpz, sco,sce,gsc,csc,vibx,viby,vibz, vibxi
     $     ,vibyi,vibzi,compress,lengthnw, pif,zif,bvx,bvy,bvz,rcb,gb
     $     ,BICM,NWICM,tavg, eqtim,lsr)


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

      DOUBLE PRECISION asr(rdim),b0r(rdim),rr_x(rdim),rr_y(rdim),
     $     rr_z(rdim), lsr(rdim)
      DOUBLE PRECISION b1c(rdim),b1s(rdim)

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
      LOGICAL tavg

      INTEGER zstart,zend,zstep

      DOUBLE PRECISION compress(rdim,adim,sdim)

      DOUBLE PRECISION trbin

      DOUBLE PRECISION sumb,sumnw
      INTEGER i,j

      DOUBLE PRECISION b2c(rdim),b2s(rdim),b3c(rdim),b3s(rdim), b4c(rdim
     $     ),b4s(rdim)

      INTEGER lengthnw(sdim)

      DOUBLE PRECISION factor,rc,rc_2p,rc_4p,zet0,wage
      INTEGER k
      INTEGER rotgeo
      
      double precision vzold(rdim,adim,sdim),vznorm

      common/pidata/pi
      save/pidata/

      common/hhflap/b2c,b2s,b3c,b3s,b4c,b4s
      save/hhflap/

      common/infnorm/vzold
      save/infnorm/
      
      common/rotdata/rotgeo
      save/rotdata/


c     Initialize total induced velocity at every blade
c     control point to zero so it isn't re-summed:
      do r0=1,nr
         do p0=1,np
            do s0=1,ns

               vbtx(r0,p0,s0)=0.D+0
               vbty(r0,p0,s0)=0.D+0
               vbtz(r0,p0,s0)=0.D+0
               vbtxin(r0,p0,s0)=0.D+0
               vbtyin(r0,p0,s0)=0.D+0
               vbtzin(r0,p0,s0)=0.D+0

            end do
         end do
      end do

      open(333,file='fort.333',status='unknown')
c     For peace of mind, initialize:

      vbx=0.D+0
      vby=0.D+0
      vbz=0.D+0

      dnp=nint(anb/(dp/DBLE(pif)))
      Nb=nint(2*pi/anb)

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

      trb=1.D+0+1.0D+0*dcy*gav/kinv
      factor=1.D+0
      k=1
c     Influenced indices:
      do r0=1,nr
         do p0=1,np
            npl=p0+np-dnp
            do s0=1,ns

c     Forcing indices:
               do r1=1,nr
                  sn=1.0D+0
                  if(rotgeo.eq.2) sn=(-1.D+0)**(r1-1)
                  do pn=p0,npl,dnp
                     p1=pn
                     if (p1 .gt. np) p1=p1-np
                     do w1=1,nw-pw !no inboard wake
                        zstart=3
                        zend=nz
                        zstep=1
                        if(nz .gt. (nint(2.D+0*2.D+0*pi/dz)*zif+1)) zend
     $                       =(nint(2.D+0*2.D+0*pi/dz)*zif+1)

                        zstart=3!+lengthnw(nw+1-w1)
                        zend=nz
                        
                        do z1=zstart,zend,zstep
                           rc=0.00855D+0*sqrt(trb*eqtim(r1,p1,w1,z1,k)
     $                          *Nb/om)

                           call vbind(r0,p0,s0,r1,p1,w1,z1,np,nz,nk,dz,
     $                          gv,rc,lin,trbin,om,rc0,cpx,cpy,cpz, pcx
     $                          ,pcy,pcz,vbx,vby,vbz,zif,Nb)

c     Sum velocity contributions for total
c     induced velocity at each control point:
                           vbtx(r0,p0,s0)=vbtx(r0,p0,s0)+vbx*sn
                           vbty(r0,p0,s0)=vbty(r0,p0,s0)+vby*sn
                           vbtz(r0,p0,s0)=vbtz(r0,p0,s0)+vbz*sn
                           vbtxin(r0,p0,s0)=vbtxin(r0,p0,s0)+vbx*sn
                           vbtyin(r0,p0,s0)=vbtyin(r0,p0,s0)+vby*sn
                           vbtzin(r0,p0,s0)=vbtzin(r0,p0,s0)+vbz*sn
                           if((s0 .eq. 15) .and. (p1 .eq. 1) .and. (p0
     $                          .eq. 1)) then
                              write(333,*) z1, vbz
                           endif
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      vbx=0.D+0
      vby=0.D+0
      vbz=0.D+0

c     Shed circulation effects:
      if (sco .eq. 'y') then
         if (nw-pw  .gt.  1 ) then !no inboard wake
c     Influenced indices:
            do r0=1,nr
               do p0=1,np
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
                              if(nz .gt. (nint(2.D+0*2.D+0*pi/dz)*zif+1)
     $                             )zend=(nint(2.D+0*2.D+0*pi/dz)*zif+1)
                              do z1=zstart,zend,zstep
                                 rc=0.00855D+0*sqrt(trb*eqtim(r1,p1,w1
     $                                ,z1,1)*Nb/om)

                                 call sbind(r0,p0,s0,r1,p1,w1,z1,np,nk
     $                                ,dz,om,rc0,gv,rc,lin,trbin,sce,gsc
     $                                ,csc,cpx,cpy,cpz,pcx,pcy,pcz,vbx
     $                                ,vby,vbz,zif)
                                 
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
         do p0=1,np
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
     $                       *compress(r1,p1,s1)
c                        vbx=vbz*sin(-asr(r1))
c                        vbz=vbz*cos(-asr(r1))
C/MR              Adding lateral shaft tilt here, though not sure
                        vbx=vbz*sin(-asr(r1))*cos(lsr(r1))
                        vbz=vbz*cos(-asr(r1))*cos(lsr(r1))                        
                        vbtxin(r0,p0,s0)=vbtxin(r0,p0,s0)+vbx
                        vbtyin(r0,p0,s0)=vbtyin(r0,p0,s0)+vby
                        vbtzin(r0,p0,s0)=vbtzin(r0,p0,s0)+vbz
                     end do
                  end do
               end do
            end do
         end do
      end do

c     Fixed to rotating frame transform:
      do r=1,nr
         sn=1.D+0
         if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r-1)
         do p=1,np
            az=DBLE(p-1)*dp/DBLE(pif)
            b0=b0r(r)+b1c(r)*cos(az)+b1s(r)*sin(az)
            az=sn*(az + ip(r))

C           t11= cos(-asr(r))*cos(az)*cos(b0)+ sin(-asr(r))*sin(b0)
C           t12= sin(az)*cos(b0)
C           t13=-sin(-asr(r))*cos(az)*cos(b0)+ cos(-asr(r))*sin(b0)
C
C           t21=-cos(-asr(r))*sin(az)*sn
C           t22= cos(az)*sn
C           t23= sin(-asr(r))*sin(az)*sn
C
C           t31=-cos(-asr(r))*cos(az)*sin(b0)+ sin(-asr(r))*cos(b0)
C           t32=-sin(az)*sin(b0)
C           t33= sin(-asr(r))*cos(az)*sin(b0)+ cos(-asr(r))*cos(b0)

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
               vibx(r,p,s)=   t11*vbtx(r,p,s)+t12*vbty(r,p,s)+t13*vbtz(r
     $              ,p,s)
               viby(r,p,s)=+( t21*vbtx(r,p,s)+t22*vbty(r,p,s)+t23*vbtz(r
     $              ,p,s) )
               vibz(r,p,s)=   t31*vbtx(r,p,s)+t32*vbty(r,p,s)+t33*vbtz(r
     $              ,p,s)
               vibxi(r,p,s)=   t11*vbtxin(r,p,s) +t12*vbtyin(r,p,s) +t13
     $              *vbtzin(r,p,s)
               vibyi(r,p,s)=+( t21*vbtxin(r,p,s) +t22*vbtyin(r,p,s) +t23
     $              *vbtzin(r,p,s) )
               vibzi(r,p,s)=   t31*vbtxin(r,p,s) +t32*vbtyin(r,p,s) +t33
     $              *vbtzin(r,p,s)
            end do

         end do
      end do

      open(71,file='VBX.dat',status='unknown')
      open(72,file='VBY.dat',status='unknown')
      open(73,file='VBZ.dat',status='unknown')

      do r=1,nr
         
         do p=1,1
            do s=ns,1,-1
               write(71,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r)) -(cpz(r
     $              ,p,s)-rr_z(r))/rad*sin(-asr(r)), vibx(r,p,s),vibxi(r
     $              ,p,s)
               write(72,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r)) -(cpz(r
     $              ,p,s)-rr_z(r))/rad*sin(-asr(r)), viby(r,p,s),vibyi(r
     $              ,p,s)
               write(73,314)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
     $              -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)), vibz(r,p,s)
     $              ,vibzi(r,p,s)
            end do
         end do

         write(71,*)
         write(72,*)
         write(73,*)

         do p=1+int(pi/(dp/DBLE(pif))),1+int(pi/(dp/DBLE(pif)))
            do s=1,ns
               write(71,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r)) -(cpz(r
     $              ,p,s)-rr_z(r))/rad*sin(-asr(r)), vibx(r,p,s),vibxi(r
     $              ,p,s)
               write(72,*)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r)) -(cpz(r
     $              ,p,s)-rr_z(r))/rad*sin(-asr(r)), viby(r,p,s),vibyi(r
     $              ,p,s)
               write(73,314)(cpx(r,p,s)-rr_x(r))/rad*cos(-asr(r))
     $              -(cpz(r,p,s)-rr_z(r))/rad*sin(-asr(r)), vibz(r,p,s)
     $              ,vibzi(r,p,s)
            end do
         end do
         
         write(71,*)
         write(72,*)
         write(73,*)

C/MR  The transformation from hub to shaft axes for x does not change with 
C/MR  the lateral shaft tilt, but y does  
C/MR  TO DO       
         do p=int(pi/(2.D+0*(dp/DBLE(pif))))+1, int(pi/(2.D+0*(dp
     $        /DBLE(pif))))+1
            do s=ns,1,-1
               write(71,*)(cpy(r,p,s)-rr_y(r))/rad, vibx(r,p,s),vibxi(r
     $              ,p,s)
               write(72,*)(cpy(r,p,s)-rr_y(r))/rad, viby(r,p,s),vibyi(r
     $              ,p,s)
               write(73,314)(cpy(r,p,s)-rr_y(r))/rad, vibz(r,p,s)
     $              ,vibzi(r,p,s)
            end do
         end do
         write(71,*)
         write(72,*)
         write(73,*)
         do p=int(pi/(2.D+0*(dp/DBLE(pif))))+1+int(pi/(dp/DBLE(pif))),
     $        int(pi/(2.D+0*(dp/DBLE(pif))))+1+int(pi/(dp/DBLE(pif)))
            do s=1,ns
               write(71,*)(cpy(r,p,s)-rr_y(r))/rad, vibx(r,p,s),vibxi(r
     $              ,p,s)
               write(72,*)(cpy(r,p,s)-rr_y(r))/rad, viby(r,p,s),vibyi(r
     $              ,p,s)
               write(73,314)(cpy(r,p,s)-rr_y(r))/rad, vibz(r,p,s)
     $              ,vibzi(r,p,s)
            end do
         end do

         write(71,*)
         write(72,*)
         write(73,*)
         write(71,*)
         write(72,*)
         write(73,*)

      end do

      close(71)
      close(72)
      close(73)

314   format(3(E12.6,1X))
      open(74,file='VBZ_PSI.dat',status='unknown')


      vznorm=0.0D0
      do r=1,nr
         do p=1,np
            do s=1,ns

               write(74,315)cpx(r,p,s),cpy(r,p,s),cpz(r,p,s), vibxi(r,p
     $              ,s),vibyi(r,p,s),vibzi(r,p,s),vibz(r,p,s)

               vznorm=vznorm+(vibzi(r,p,s)-vzold(r,p,s))**2
               vzold(r,p,s)=vibzi(r,p,s)
            end do
            write(74,*)
            write(74,*)
         end do
         write(74,*)
      end do
      
      write(*,*) "Inflow norm", sqrt(vznorm/(DBLE(nr*np*ns)))
      close(74)
 315  format(7(E12.6,1X))

 10   format(30(f10.5,1x))

      close(333)
      return
      end 
