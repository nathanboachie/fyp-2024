c----------------------------------------------------------------------|
c      Copyright S
c
c----------------------------------------------------------------------|
c
      subroutine strainab(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     &                    gv,rc0,kinv,dcy,om,gav,
     &                    ptx,pty,ptz, eqtim)

      IMPLICIT none
c     
      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'
c     
      INTEGER nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,
     &     r,p,w,z,n,s,k,pm1
      
      DOUBLE PRECISION ip(rdim), dp, dz, kinv, dcy,om, gav

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim),
     &       rc0(rdim,adim,wdim,kdim),
     &       ptx(rdim,adim,wdim,zdim),
     &       pty(rdim,adim,wdim,zdim),
     &       ptz(rdim,adim,wdim,zdim),
     &       eqtim(rdim,adim,wdim,zdim,kdim)

      DOUBLE PRECISION trb,wage,zet0,factor,rc

      DOUBLE PRECISION v1x,v2x,v1y,v2y,v1z,v2z,vec1,vec2,rat

      factor=1.0D+0
      do r=1,nr
         do p=1,np
            pm1=p-1
            if (pm1 .le. 0) pm1=pm1+np
            do w=1,nw
               do z=3,nzt
                  v1x=ptx(r,pm1,w,z-1)-ptx(r,pm1,w,z-2)
                  v1y=pty(r,pm1,w,z-1)-pty(r,pm1,w,z-2)
                  v1z=ptz(r,pm1,w,z-1)-ptz(r,pm1,w,z-2)
                  v2x=ptx(r,p,w,z)-ptx(r,p,w,z-1)
                  v2y=pty(r,p,w,z)-pty(r,p,w,z-1)
                  v2z=ptz(r,p,w,z)-ptz(r,p,w,z-1)
                  vec1=(v1x*v1x+v1y*v1y+v1z*v1z)**0.5D+0
                  vec2=(v2x*v2x+v2y*v2y+v2z*v2z)**0.5D+0
                  rat=(vec1/vec2)**0.5D+0
                  do k=1,nk
                     wage=DBLE(z-1)*dz/DBLE(zif)
                     trb=1.D+0+dcy*factor*abs(gv(r,pm1,w,z-1,1)/kinv)
                     trb=1.D+0+dcy*factor*gav/kinv
                     rc=0.00855D+0*rat*
     &                sqrt(trb*eqtim(r,pm1,w,z-1,k)*Nb/om)
                     zet0=(rc/0.00855D+0)**2*om/trb/Nb
                     eqtim(r,p,w,z,k)=zet0+dz
                  end do
               end do
            end do
         end do
      end do

      return
      end
