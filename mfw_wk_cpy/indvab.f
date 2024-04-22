C-$	$Id: indvab.f,v 1.3.2.3 2005/06/11 15:57:37 shreyas Exp $	

c    <Copyright (c) 1997 by Ashish Bagai>


      subroutine indvab(nr,np,nw,pw,nz,nzt,nk,ns,pif,zif,dp,dz,rad, om
     $     ,anb,lin,dcy,gav,kinv,gv,rc0,sco,sce,gsc,csc, rcb,gb,ip,b0r
     $     ,asr,rr_x,rr_y,rr_z, t0r,t1c,t1s,b1c,b1s, spx,spy,spz,px,py
     $     ,pz,bvx,bvy,bvz, vitx,vity,vitz,cpx,cpy,cpz, eqtim)


c     Induced velocity calculations -- free-vortex/free-vortex and
c     bound-vortex/free-vortex interactions for total induced velocity
c     at free collocation points in wake.  Pass geometries, initial core
c     radii, strengths in.  Total induced velocities returned.

c     r => rotor index
c     p => azimuthal index
c     w => vortex filament index
c     z => vortex element index
c     s => blade segment index


      IMPLICIT none


      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER  nr,np,nw,pw,nz,nzt,nk,ns,pif,zif,pn,
     &     r0,p0,w0,z0,r1,p1,w1,z1,s1,npl,dnp

      DOUBLE PRECISION dp,dz,om,anb,lin,dcy,gav,kinv,
     &     vix,viy,viz,trb,ip(rdim),rcb,sce,wage

      DOUBLE PRECISION b0r(rdim),asr(rdim),
     &     spx(sdim),spy(sdim),spz(sdim),
     &     rr_x(rdim),rr_y(rdim),rr_z(rdim)

      DOUBLE PRECISION b1c(rdim),b1s(rdim)

      DOUBLE PRECISION bvx(rdim,adim,sdim), bvy(rdim,adim,sdim),
     $     bvz(rdim,adim,sdim), gb(rdim,adim,sdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), rc0(rdim,adim,wdim
     $     ,kdim), px(rdim,adim,wdim,zdim), py(rdim,adim,wdim,zdim),
     $     pz(rdim,adim,wdim,zdim), vitx(rdim,adim,wdim,zdim), vity(rdim
     $     ,adim,wdim,zdim), vitz(rdim,adim,wdim,zdim), eqtim(rdim,adim
     $     ,wdim,zdim,kdim)

      CHARACTER*1 sco

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim),
     &     csc(rdim,adim,wdim,zdim,kdim)

      INTEGER lr,lp,lk

      INTEGER itmp,PsiO,Nb

      DOUBLE PRECISION cpx(rdim,adim,sdim), cpy(rdim,adim,sdim),
     $     cpz(rdim,adim,sdim)

      DOUBLE PRECISION rad,ratio,pi,aza,theta,viz1,viz2

      INTEGER s0

      LOGICAL bound,deltam,noindv

      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim)

      DOUBLE PRECISION factor,rc,rc_2p,rc_4p,zet0
      INTEGER k

      DOUBLE PRECISION rcore(rdim,adim,wdim,zdim) ,squire(rdim,adim,wdim
     $     ,zdim) ,osquire(rdim,adim,wdim,zdim)

      INTEGER pm1

      DOUBLE PRECISION sn
      double precision nx1(sdim),ny1(sdim),nz1(sdim),nwx(rdim,adim,sdim)
     $     ,nwy(rdim,adim,sdim),nwz(rdim,adim,sdim)



      common/coredata/rcore,squire,osquire
      save/coredata/

      common/nearwake/nx1,ny1,nz1
      save/nearwake/


      pi=4.D+0*atan(1.D+0)


c     Initialize total induced velocity at every wake
c     collocation point to zero so it isn't re-summed:
      do r0=1,nr
         do p0=1,np
            do w0=1,nw
               do z0=1,nz
                  
                  vitx(r0,p0,w0,z0)=0.D+0
                  vity(r0,p0,w0,z0)=0.D+0
                  vitz(r0,p0,w0,z0)=0.D+0

               end do
            end do
         end do
      end do


      dnp=nint(anb/(dp/DBLE(pif)))

      Nb=nint(2*pi/anb)

c     Rotating to fixed frame coordinate transformation:

      call rtof(nr,np,ns,pif,dp,ip,b0r,b1c,b1s,asr, rr_x,rr_y,rr_z,spx
     $     ,spy,spz, bvx,bvy,bvz)
      call rtof(nr,np,ns,pif,dp,ip,b0r,b1c,b1s,asr, rr_x,rr_y,rr_z,nx1
     $     ,ny1,nz1, nwx,nwy,nwz)


c     For peace of mind, initialize:
      vix=0.D+0
      viy=0.D+0
      viz=0.D+0

      deltam=.false.
      if(deltam) then
c     Influenced indices:
         do r0=1,nr
            do p0=1,np,pif
               npl=p0+np-dnp

               do w0=1,nw !  calculate the "prescribed" inboard sheet 
                  do z0=1,nz,zif !note 1 instead of 1+zif!
                     PsiO=p0-(z0-1)
                     do itmp=1,nzt
                        if(PsiO.le.0) PsiO=PsiO+np
                     end do
                     aza=ip(r0)+DBLE(p0-1)*dp/DBLE(pif) - ip(r0)
                     s0=(ns+1-w0)

                     do r1=r0,r0
                        sn=(-1.D+0)**(r1-1)
                        do p1=PsiO,PsiO,1
                           do s1=1,ns !s0,s0

                              aza=ip(r1)+DBLE(PsiO-1)*dp/DBLE(pif) -
     $                             ip(r1)
                              call momentum(r0,PsiO,s0,ns,r1,p1,s1,rcb
     $                             ,gb,rad,aza,t0r,t1c,t1s,cpx,cpy,cpz
     $                             ,bvx,bvy,bvz,vix,viy,viz)

c     Sum velocity contributions for total induced
c     velocity at each collocation point:
                              vitx(r0,p0,w0,z0)=vitx(r0,p0,w0,z0)+vix*sn
                              vity(r0,p0,w0,z0)=vity(r0,p0,w0,z0)+viy*sn
                              vitz(r0,p0,w0,z0)=vitz(r0,p0,w0,z0)+viz*sn

                           end do
                        end do
                     end do

                  end do
               end do
            end do
         end do
      endif                     !to compute momentum addition


c     For peace of mind, initialize:
      vix=0.D+0
      viy=0.D+0
      viz=0.D+0


c     Free-vortex/free-vortex interactions:

c     Trailed circulation effects:


c     Influenced indices:
      do r0=1,nr
         do p0=1,np,pif
            npl=p0+np-dnp
            do w0=1,nw-pw
               do z0=1,nz,zif

c     Forcing indices:
                  do r1=1,nr
                     sn=(-1.D+0)**(r1-1)
                     do pn=p0,npl,dnp
                        p1=pn
                        if (p1 .gt. np) p1=p1-np
                        do w1=1,nw-pw

                           gav=0.D+0
                           do lr=1,nr
                              do lp=1,np
                                 do lk=1,nk
                                    gav=gav+abs(gv(lr,lp,w1,1,lk))
                                 end do
                              end do
                           end do
                           gav=gav/DBLE(nr*np)
                           trb=1.D+0+dcy*gav/kinv
                           k=1
                           do z1=(1+zif),nzt,zif

                              wage=DBLE(z1-1)*dz/DBLE(zif)
                              trb=1.D+0+dcy*osquire(r1,p1,w1,z1)
     $                             *abs(gv(r1,p1,w1,z1,1)/kinv)
                              trb=1.D+0+dcy*osquire(r1,p1,w1,z1)*gav
     $                             /kinv
                              trb=1.D+0+dcy*factor*gav/kinv
                              pm1=p1-1
                              if(pm1.le.0)pm1=pm1+np
                              zet0=(rcore(r1,pm1,w1,z1-1)/.00855D+0)**2
     $                             *om/trb/Nb-(wage-dz)

                              rc=0.00855D+0*sqrt(trb*(wage+zet0)*Nb/om)
                              factor=1.D+0
c-    S Modified core radius approach...
                              rc=0.00855D+0*sqrt(trb*eqtim(r1,p1,w1,z1,k
     $                             )*Nb/om)

                              call vvind(r0,p0,w0,z0,r1,p1,w1,z1,np,nz
     $                             ,nzt,nk,pif,zif,dz,om,anb,gv,rc,rc0
     $                             ,lin,trb,px,py,pz,vix,viy,viz,Nb
     $                             ,factor)

                              rcore(r1,p1,w1,z1)=rc
                              squire(r1,p1,w1,z1)=factor
                              
c     Sum velocity contributions for total induced
c     velocity at each collocation point:
                              vitx(r0,p0,w0,z0)=vitx(r0,p0,w0,z0)+vix*sn
                              vity(r0,p0,w0,z0)=vity(r0,p0,w0,z0)+viy*sn
                              vitz(r0,p0,w0,z0)=vitz(r0,p0,w0,z0)+viz*sn

                           end do
                        end do
                     end do
                  end do

               end do
            end do
         end do
      end do

c     For peace of mind, reinitialize:
      vix=0.D+0
      viy=0.D+0
      viz=0.D+0


c     Shed circulation effects:

      if (sco .eq. 'y') then
         if (nw-pw  .gt.  1 ) then

c     Influenced indices:
            do r0=1,nr
               do p0=1,np,pif
                  npl=p0+np-dnp
                  do w0=1,nw-pw
                     do z0=1,nz,zif

c     Forcing indices:
                        do r1=1,nr
                           do pn=p0,npl,dnp
                              p1=pn
                              if (p1 .gt. np) p1=p1-np
                              do w1=2,nw-pw
                                 do z1=(1+zif),nzt-zif,zif 

                                    factor=1.D+0
                                    k=1
                                    trb=1.D+0+dcy*factor*gav/kinv
                                    rc=0.00855D+0*sqrt(trb*eqtim(r1,p1
     $                                   ,w1,z1,k)*Nb/om)

                                    call svind(r0,p0,w0,z0,r1,p1,w1,z1
     $                                   ,np,nk,pif,zif,dz,om,anb,gv,rc
     $                                   ,rc0,lin,trb,sce,gsc,csc,px,py
     $                                   ,pz,vix,viy,viz)

c     Sum velocity contributions for total induced
c     velocity at each collocation point:
                                    vitx(r0,p0,w0,z0)=vitx(r0,p0,w0,z0)
     $                                   +vix
                                    vity(r0,p0,w0,z0)=vity(r0,p0,w0,z0)
     $                                   +viy
                                    vitz(r0,p0,w0,z0)=vitz(r0,p0,w0,z0)
     $                                   +viz

                                 end do
                              end do
                           end do
                        end do

                     end do
                  end do
               end do
            end do

         endif
      endif


c     For peace of mind, reinitialize:
      vix=0.D+0
      viy=0.D+0
      viz=0.D+0
      bound=.true.
      if(bound) then


c     Bound-vortex/free-vortex interactions:

c     Rotating to fixed frame coordinate transformation:
         call rtof(nr,np,ns,pif,dp,ip,b0r,b1c,b1s,asr, rr_x,rr_y,rr_z
     $        ,spx,spy,spz, bvx,bvy,bvz)

c     Influenced indices:
         do r0=1,nr
            do p0=1,np,pif
               npl=p0+np-dnp
               do w0=1,nw-pw
                  do z0=1+zif,nz,zif !note 1 instead of 1+zif!

c     Forcing indices:
                     do r1=1,nr
                        sn=(-1.D+0)**(r1-1)
                        do pn=p0,p0!npl,dnp
                           p1=pn
                           if (p1 .gt. np) p1=p1-np
                           do s1=1,ns

                              call bvind(r0,p0,w0,z0,r1,p1,s1,rcb,gb, px
     $                             ,py,pz,bvx,bvy,bvz,vix,viy,viz)

                              vitx(r0,p0,w0,z0)=vitx(r0,p0,w0,z0)+vix*sn
                              vity(r0,p0,w0,z0)=vity(r0,p0,w0,z0)+viy*sn
                              vitz(r0,p0,w0,z0)=vitz(r0,p0,w0,z0)+viz*sn

                              call nwind(r0,p0,w0,1+zif,r1,p1,s1,rcb,gb
     $                             ,px,py,pz,bvx,bvy,bvz,nwx,nwy,nwz,vix
     $                             ,viy,viz)

                              vitx(r0,p0,w0,z0)=vitx(r0,p0,w0,z0)+vix*sn
                              vity(r0,p0,w0,z0)=vity(r0,p0,w0,z0)+viy*sn
                              vitz(r0,p0,w0,z0)=vitz(r0,p0,w0,z0)+viz*sn


                           end do
                        end do
                     end do

                  end do
               end do
            end do
         end do


      endif                     !bound

c     Induced velocity interpolation for ersatz-free points:

      call vinterp(nr,np,nw,pw,nz,pif,zif,vitx,vity,vitz)


 1014 continue


      do r0=1,nr
         do w0=1,nw
            do z0=1,nzt
               do p0=1,np
                  pm1=p0-1
                  if(pm1.le.0)pm1=pm1+np
                  osquire(r0,p0,w0,z0)=squire(r0,p0,w0,z0)
                  if(z0.gt.1) then
                     if(osquire(r0,pm1,w0,z0-1).gt.osquire(r0,p0,w0,z0))
     $                    osquire(r0,p0,w0,z0)=osquire(r0,pm1,w0,z0-1)
                  endif

                  if(z0.gt.nzt) then
                     osquire(r0,p0,w0,z0)=osquire(r0,pm1,w0,z0-1)
                  endif

               end do
            end do
         end do
      end do



      return
      end
