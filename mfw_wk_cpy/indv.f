C-$	$Id: indv.f,v 1.2.4.2 2005/06/11 15:57:37 shreyas Exp $	

      subroutine indv(nr,np,nw,pw,nz,nzt,nk,ns,pif,zif,dp,dz,rad, om,anb
     $     ,lin,dcy,gav,kinv,gv,rc0,sco,sce,gsc,csc, rcb,gb,ip,b0r,asr
     $     ,lsr,rr_x,rr_y,rr_z, t0r,t1c,t1s,beta, spx,spy,spz,px,py,pz
     $     ,bvx,bvy,bvz, vitx,vity,vitz,p0,cpx,cpy,cpz, eqtim)

C-$   Compute induced velocities at wake control points for PC2B 


      IMPLICIT none

      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER  nr,np,nw,pw,nz,nzt,nk,ns,pif,zif,pn, r0,p0,w0,z0,r1,p1,w1
     $     ,z1,s1,npl,dnp

      DOUBLE PRECISION dp,dz,om,anb,lin,dcy,gav,kinv, vix,viy,viz,trb
     $     ,ip(rdim),rcb,sce,wage

      DOUBLE PRECISION b0r(rdim),asr(rdim), spx(sdim),spy(sdim),spz(sdim
     $     ), rr_x(rdim),rr_y(rdim),rr_z(rdim),lsr(rdim)

      DOUBLE PRECISION bvx(rdim,adim,sdim), bvy(rdim,adim,sdim),
     $     bvz(rdim,adim,sdim), gb(rdim,adim,sdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), rc0(rdim,adim,wdim
     $     ,kdim), px(rdim,adim,wdim,zdim), py(rdim,adim,wdim,zdim),
     $     pz(rdim,adim,wdim,zdim), vitx(rdim,adim,wdim,zdim), vity(rdim
     $     ,adim,wdim,zdim), vitz(rdim,adim,wdim,zdim), eqtim(rdim,adim
     $     ,wdim,zdim,kdim)

      CHARACTER*1 sco

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim), csc(rdim,adim,wdim
     $     ,zdim,kdim)

      INTEGER lr,lp,lk

      INTEGER itmp,PsiO,Nb

      DOUBLE PRECISION cpx(rdim,adim,sdim), cpy(rdim,adim,sdim),
     $     cpz(rdim,adim,sdim)

      DOUBLE PRECISION rad,ratio,pi,aza,theta

      INTEGER s0

      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim)
      DOUBLE PRECISION beta(rdim,adim)

      DOUBLE PRECISION factor,rc,rc_2p,rc_4p,zet0
      INTEGER k

      DOUBLE PRECISION sn
      INTEGER rotgeo

      double precision nx1(sdim),ny1(sdim),nz1(sdim),nwx(rdim,adim,sdim)
     $     ,nwy(rdim,adim,sdim),nwz(rdim,adim,sdim)

      LOGICAL bound, deltam

      double precision vcx_ml(rdim*adim*zdim*wdim)
      double precision vcy_ml(rdim*adim*zdim*wdim)
      double precision vcz_ml(rdim*adim*zdim*wdim)
      INTEGER idx

      common/nearwake/nx1,ny1,nz1
      save/nearwake/

      common/rotdata/rotgeo
      save/rotdata/

      common/indv_ml/vcx_ml,vcy_ml,vcz_ml
      save/indv_ml/

      pi=4.D+0*atan(1.D+0)

c >>  Initialize total induced velocity at every wake  <<
c >>  collocation point to zero so it isn't re-summed: <<
      do r0=1,nr
         do w0=1,nw
            do z0=1,nz
               
               vitx(r0,p0,w0,z0)=0.D+0
               vity(r0,p0,w0,z0)=0.D+0
               vitz(r0,p0,w0,z0)=0.D+0

            end do
         end do
      end do

      if(.false.) then
      do r0=1,nr
         npl=p0+np-dnp
         do w0=1,nw         
            do z0=1,nz,zif   
               PsiO=p0-(z0-1)
               do itmp=1,nzt
                  if(PsiO.le.0) PsiO=PsiO+np
               end do
               aza=ip(r0)+DBLE(p0-1)*dp/DBLE(pif) - ip(r0)
               s0=(ns+1-w0)
               do r1=r0,r0
                  sn=1.D+0
                 if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r1-1)
                  do p1=PsiO,PsiO,1
                     do s1=1,ns

                        aza=ip(r1)+DBLE(PsiO-1)*dp/DBLE(pif) - ip(r1)
                        call momentum(r0,PsiO,s0,ns,r1,p1,s1,rcb,gb,rad
     $                       ,aza,t0r,t1c,t1s,cpx,cpy,cpz,bvx,bvy,bvz
     $                       ,vix,viy,viz)
                        vitx(r0,p0,w0,z0)=vitx(r0,p0,w0,z0)+vix*sn
                        vity(r0,p0,w0,z0)=vity(r0,p0,w0,z0)+viy*sn
                        vitz(r0,p0,w0,z0)=vitz(r0,p0,w0,z0)+viz*sn
                     end do
                  end do
               end do

            end do
         end do
      end do
      endif

c     For peace of mind, initialize:

      vix=0.D+0
      viy=0.D+0
      viz=0.D+0


      dnp=nint(anb/(dp/DBLE(pif)))

      Nb=nint(2*pi/anb)

c >>  Free-vortex/free-vortex interactions:  <<
c >>  Trailed circulation effects:  <<

c >>  Influenced indices:   <<
      do r0=1,nr
         npl=p0+np-dnp
         do w0=1,nw-pw
            do z0=1,nz,zif

c     >>     Forcing indices:  <<
               do r1=1,nr
                  sn=1.D+0
                  if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r1-1)
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
                        do z1=(1+zif),nzt,zif
                           k=1
                           rc=0.00855D+0*sqrt(trb*eqtim(r1,p1,w1,z1,k)
     $                          *Nb/om)
                           call vvind(r0,p0,w0,z0,r1,p1,w1,z1,np,nz,nzt
     $                          ,nk,pif,zif,dz,om,anb,gv,rc,rc0,lin,trb
     $                          ,px,py,pz,vix,viy,viz,Nb,factor)

c >>        Sum velocity contributions for total induced <<
c >>        velocity at each collocation point:          <<
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

c >>  For peace of mind, reinitialize: <<

      vix=0.D+0
      viy=0.D+0
      viz=0.D+0

c >>  Shed circulation effects: <<

      if (sco .eq. 'y') then
         if (nw-pw  .gt.  1 ) then

c     >>  Influenced indices: <<
            do r0=1,nr
               npl=p0+np-dnp
               do w0=1,nw-pw
                  do z0=1,nz,zif

c     >>     Forcing indices: <<
                     do r1=1,nr
                        do pn=p0,npl,dnp
                           p1=pn
                           if (p1 .gt. np) p1=p1-np
                           do w1=2,nw-pw
                              do z1=(1+zif),nzt-zif,zif
                                 factor=1.D+0
                                 k=1
                                 trb=1.D+0+dcy*factor*gav/kinv
                                 rc=0.00855D+0*sqrt(trb*eqtim(r1,p1,w1
     $                                ,z1,k)*Nb/om)


                                 call svind(r0,p0,w0,z0,r1,p1,w1,z1,np
     $                                ,nk,pif,zif,dz,om,anb,gv,rc,rc0
     $                                ,lin,trb,sce,gsc,csc,px,py,pz,vix
     $                                ,viy,viz)

c >>        Sum velocity contributions for total induced <<
c >>        velocity at each collocation point:          <<
                                 vitx(r0,p0,w0,z0)=vitx(r0,p0,w0,z0)+vix
                                 vity(r0,p0,w0,z0)=vity(r0,p0,w0,z0)+viy
                                 vitz(r0,p0,w0,z0)=vitz(r0,p0,w0,z0)+viz

                              end do
                           end do
                        end do
                     end do

                  end do
               end do
            end do

         endif
      endif

c >>  For peace of mind, reinitialize: <<

      vix=0.D+0
      viy=0.D+0
      viz=0.D+0

c >>  Bound-vortex/free-vortex interactions: <<

c >>  Influenced indices: <<
      if (.true.) then
      do r0=1,nr
         npl=p0+np-dnp
         do w0=1,nw-pw
            do z0=1+zif,nz,zif

c     >>     Forcing indices: <<
               do r1=1,nr
                  sn=1.D+0
                  if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r1-1)
                  do pn=p0,p0!npl,dnp
                     p1=pn
                     if (p1 .gt. np) p1=p1-np
                     call rtofpsi(nr,np,ns,pif,dp,ip,beta,asr,lsr,rr_x
     $                    ,rr_y,rr_z,spx,spy,spz, bvx,bvy,bvz,p1)
                     call rtofpsi(nr,np,ns,pif,dp,ip,beta,asr,lsr,rr_x
     $                    ,rr_y,rr_z,nx1,ny1,nz1,nwx,nwy,nwz,p1)

                     do s1=1,ns

                        call bvind(r0,p0,w0,z0,r1,p1,s1,rcb,gb, px,py,pz
     $                       ,bvx,bvy,bvz,vix,viy,viz)

c >>        Sum velocity contributions for total induced <<
c >>        velocity at each collocation point:          <<
                        vitx(r0,p0,w0,z0)=vitx(r0,p0,w0,z0)+vix*sn
                        vity(r0,p0,w0,z0)=vity(r0,p0,w0,z0)+viy*sn
                        vitz(r0,p0,w0,z0)=vitz(r0,p0,w0,z0)+viz*sn

                        call nwind(r0,p0,w0,1+zif,r1,p1,s1,rcb,gb,px,py
     $                       ,pz,bvx,bvy,bvz,nwx,nwy,nwz,vix,viy,viz)

                        vitx(r0,p0,w0,z0)=vitx(r0,p0,w0,z0)+vix*sn
                        vity(r0,p0,w0,z0)=vity(r0,p0,w0,z0)+viy*sn
                        vitz(r0,p0,w0,z0)=vitz(r0,p0,w0,z0)+viz*sn

                     end do
                  end do
               end do

            end do
         end do
      end do
      endif
      return

      idx = 1
c >>  Influenced indices, Fuselage velocities:   <<
      do r0=1,nr
         do w0=1,nw
            do z0=1,nz
            vitx(r0,p0,w0,z0)=vitx(r0,p0,w0,z0)+vcx_ml(idx)
            vity(r0,p0,w0,z0)=vity(r0,p0,w0,z0)+vcy_ml(idx)
            vitz(r0,p0,w0,z0)=vitz(r0,p0,w0,z0)+vcz_ml(idx)
            idx=idx+1
            end do
         end do
      end do
      end
