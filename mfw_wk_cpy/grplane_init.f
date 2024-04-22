
      subroutine grplane_init(nz,rad,om,dz)

      implicit none

      include 'grdim.inc'

      integer i,j,nz

      double precision rad,dx,dy,dz,om
      double precision xmin,ymin,xmax,ymax
      double precision grxco(grxdim), gryco(grydim),grvz(grxdim,grydim)
     $     ,grmatrix(grxdim*grydim,grxdim*grydim)

      common/grplanedata/grxco,gryco,grvz
      save/grplanedata/

      xmin=-5.0D0*rad
      xmax=(5.0D0+mu*om*nz*dp)*rad
      ymin=-5.0D0*rad
      ymax=5.0D0*rad

      dx=(xmax-xmin)/(DBLE(grxdim-1))
      dy=(ymax-ymin)/(DBLE(grydim-1))
      
  
       grxco(1)=xmin
       do i=2,grxdim
          grxco(i)=grxco(i-1)+dx
       enddo
      
      gryco(1)=ymin
      do i=2,grydim
         gryco(i)=gryco(i-1)+dy
      enddo

      do j=1,grydim
         do i=1,grxdim
            grvz(i,j)=0.0D0
         enddo
      enddo

      return
      end


