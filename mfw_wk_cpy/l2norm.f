c----------------------------------------------------------------------|
c    <Copyright (c) 1997 by Ashish Bagai>
c     Copyright violation -MrB changed MRS to RMS!
c----------------------------------------------------------------------|
c
      subroutine l2norm(nr,np,nw,nz,pox,poy,poz,pcx,pcy,pcz,
     &                  gs,gsv,gsno,rms)
c     Wake convergence RMS calculations.
c----------------------------------------------------------------------|
c-mb check this thing- corrected for that
c-mb i think R M S means root of mean square
c-mb here it's more like MRS
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'rdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'
c
      INTEGER nr,np,nw,nz,gs,gsno,r,p,w,z,nzforrms
c
      DOUBLE PRECISION dx,dy,dz
c
      DOUBLE PRECISION gsv(gdim,gdim),rms(rdim,wdim)
c
      DOUBLE PRECISION pcx(rdim,adim,wdim,zdim),
     &       pcy(rdim,adim,wdim,zdim),
     &       pcz(rdim,adim,wdim,zdim),
     &       pox(rdim,adim,wdim,zdim),
     &       poy(rdim,adim,wdim,zdim),
     &       poz(rdim,adim,wdim,zdim)
c
      DOUBLE PRECISION maxrms,error
      INTEGER maxr,maxp,maxw,maxz
c----------------------------------------------------------------------|
c
      nzforrms=nz
c
      do r=1,nr
      do w=1,nw
         rms(r,w)=0.D+0
      end do
      end do
c
      maxrms=0.D+0
      do r=1,nr
      do p=1,np
      do w=1,nw
      do z=2,nzforrms,1  !not include the blade attachment point (:-)
         dx=pcx(r,p,w,z)-pox(r,p,w,z)
         dy=pcy(r,p,w,z)-poy(r,p,w,z)
         dz=pcz(r,p,w,z)-poz(r,p,w,z)
         error=dx**2+dy**2+dz**2
         rms(r,w)=rms(r,w)+error
         if(error.gt.maxrms) then
         maxrms=error
         maxr=r
         maxp=p
         maxw=w
         maxz=z
         endif
      end do
      end do
      end do
      end do
c----------------------------------------------------------------------|
c
      if (gs .eq. 1) then
      do r=1,nr
      do w=1,nw
c-mb this looks like MRS!
c-mb         rms(r,w)=(1.D+0/( (DBLE(np)/DBLE(gsv(gsno-1,3)))*
c-mb     &            (DBLE(nz-1)/DBLE(gsv(gsno-1,3))+1.D+0) ))*
c-mb     &            sqrt(rms(r,w))
c-mb this is more like RMS (:-)
         rms(r,w)=sqrt( (1.D+0/( (DBLE(np)/DBLE(gsv(gsno-1,3)))*
     &                (DBLE(nzforrms-1)/DBLE(gsv(gsno-1,3))+1.D+0) ))*
     &                   rms(r,w) )
      end do
      end do
      endif
c----------------------------------------------------------------------|
c
      if (gs .ne. 1) then
      do r=1,nr
      do w=1,nw
c-mb         rms(r,w)=(1.D+0/DBLE(np*nz))*sqrt(rms(r,w))
         rms(r,w)=sqrt( (1.D+0/DBLE(np*nzforrms))*rms(r,w) )
      end do
      end do
      endif
c----------------------------------------------------------------------|
c
      write(*,'(A12,4(1X,I4))')'Max change: ',maxr,maxp,maxw,maxz
      return
      end
