C-$	$Id: strainmb.f,v 1.3.2.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine strainmb(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     $     gv,rc0,kinv,dcy,om,gav, ppx,ppy,ppz, pcx,pcy,pcz,p,prco,
     $     eqtim)

C-$   Vortex filament strain model

      IMPLICIT none

      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif, r,p,w,z,n,s,k,pm1
      
      DOUBLE PRECISION ip(rdim), dp, dz, kinv, dcy,om, gav

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), rc0(rdim,adim,wdim
     $     ,kdim), ppx(rdim,adim,wdim,zdim), ppy(rdim,adim,wdim,zdim),
     $     ppz(rdim,adim,wdim,zdim), pcx(rdim,adim,wdim,zdim), pcy(rdim
     $     ,adim,wdim,zdim), pcz(rdim,adim,wdim,zdim), eqtim(rdim,adim
     $     ,wdim,zdim,kdim)

      DOUBLE PRECISION trb,wage,zet0,factor,rc

      DOUBLE PRECISION v1x,v2x,v1y,v2y,v1z,v2z,vec1,vec2,rat

      CHARACTER*1 prco
      factor=1.0D+0
      do r=1,nr
         pm1=p-1
         if (pm1 .le. 0) pm1=pm1+np
         do w=1,nw
            do z=3,nzt
               v1x=pcx(r,pm1,w,z-1)-pcx(r,pm1,w,z-2)
               v1y=pcy(r,pm1,w,z-1)-pcy(r,pm1,w,z-2)
               v1z=pcz(r,pm1,w,z-1)-pcz(r,pm1,w,z-2)
               if(prco .eq. 'p') then
                  v2x=ppx(r,p,w,z)-ppx(r,p,w,z-1)
                  v2y=ppy(r,p,w,z)-ppy(r,p,w,z-1)
                  v2z=ppz(r,p,w,z)-ppz(r,p,w,z-1)
               elseif (prco .eq. 'c') then
                  v2x=pcx(r,p,w,z)-pcx(r,p,w,z-1)
                  v2y=pcy(r,p,w,z)-pcy(r,p,w,z-1)
                  v2z=pcz(r,p,w,z)-pcz(r,p,w,z-1)
               else
                  v2x=v1x
                  v2y=v1y
                  v2z=v1z
               endif
               vec1=(v1x*v1x+v1y*v1y+v1z*v1z)**0.5D+0
               vec2=(v2x*v2x+v2y*v2y+v2z*v2z)**0.5D+0
               rat=(vec1/vec2)**0.5D+0
               do k=1,nk
                  wage=DBLE(z-1)*dz/DBLE(zif)
                  trb=1.D+0+dcy*factor*abs(gv(r,pm1,w,z-1,1)/kinv)
                  trb=1.D+0+dcy*factor*gav/kinv
                  rc=0.00855D+0*rat*
     &             sqrt(trb*eqtim(r,pm1,w,z-1,k)*Nb/om)
                  zet0=(rc/0.00855D+0)**2*om/trb/Nb
                  eqtim(r,p,w,z,k)=zet0+dp
               end do
            end do
         end do
      end do
      return
      end
