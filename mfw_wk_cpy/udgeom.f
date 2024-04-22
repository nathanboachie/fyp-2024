C-$	$Id: udgeom.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine udgeom(nr,np,nw,nz,nzt,pox,poy,poz,ppx,ppy,ppz, pcx,pcy
     $     ,pcz,Iteratn,rbc,zbc,avgzbc,dz,zif)

C-$   Update wake geometry positions at the end of each iteration

      IMPLICIT none


      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'


      INTEGER Iteratn

      INTEGER nr,np,nw,nz,nzt,r,p,w,z

      DOUBLE PRECISION ppx(rdim,adim,wdim,zdim), ppy(rdim,adim,wdim,zdim
     $     ), ppz(rdim,adim,wdim,zdim), pcx(rdim,adim,wdim,zdim),
     $     pcy(rdim,adim,wdim,zdim), pcz(rdim,adim,wdim,zdim), pox(rdim
     $     ,adim,wdim,zdim), poy(rdim,adim,wdim,zdim), poz(rdim,adim
     $     ,wdim,zdim)


      INTEGER pp,nz1,zif
      DOUBLE PRECISION rbc(rdim,wdim),zbc(rdim,wdim),avgzbc(rdim,wdim)
      DOUBLE PRECISION sum1,sum2,sum3,sumrad,radmin,pi,dz

c     Previous iteration geometry required for L2Norm. For first
c     iteration, these were identical to initial prescribed wake.

      do r=1,nr
         do p=1,np
            do w=1,nw
               do z=1,nzt
                  if(Iteratn.eq.1 .OR. z.ne.1) then
                     pox(r,p,w,z)=pcx(r,p,w,z)
                     poy(r,p,w,z)=pcy(r,p,w,z)
                     poz(r,p,w,z)=pcz(r,p,w,z)
                  endif
                  ppx(r,p,w,z)=pcx(r,p,w,z)
                  ppy(r,p,w,z)=pcy(r,p,w,z)
                  ppz(r,p,w,z)=pcz(r,p,w,z)

               end do
            end do
         end do
      end do

      pi=4.D+0*atan(1.D+0)
      nz1=nint(2.D+0*pi/dz)*zif         !number of pts in one turn

      do r=1,nr
         do w=1,nw    

            sum1=0.D+0
            sum2=0.D+0
            sum3=0.D+0
            sumrad=0.D+0
            radmin=1.D+9

            do pp=1,np
               sum1=sum1+poz(r,pp,w,nz)
               sum2=sum2+poz(r,pp,w,nz-nz1)
               sum3=sum3+poz(r,pp,w,nz-nz1-nz1)
               sumrad=sumrad+sqrt(pox(r,pp,w,nz)**2+poy(r,pp,w,nz)**2)
               if( radmin.gt.sqrt(pox(r,pp,w,nz)**2+poy(r,pp,w,nz)**2)
     $              )radmin = sqrt(pox(r,pp,w,nz)**2+poy(r,pp,w,nz)**2) 
            end do
            zbc(r,w)=((sum1-sum2)/1.D+0)/DBLE(np)/DBLE(nz1)
            sumrad=sumrad/DBLE(np)
            rbc(r,w)=sumrad
            avgzbc(r,w)=sum1
         end do
      end do


      return
      end
