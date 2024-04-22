C-$	$Id: blade.f,v 1.4.2.2 2005/06/11 15:57:36 shreyas Exp $	

      subroutine blade(nr,ns,crd,rcb,bcx1,bcy1,spx,spy,spz,rad,taperst
     $     ,taper,dp,lengthnw,bicm,nwicm,tw,twt,t0r)

      IMPLICIT none

      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER r,nr,s,ns,i,j

      DOUBLE PRECISION rcn,crd,dpn,pi,mr1,mr2,mr3, r1x,r1y ,r1z,r2x,r2y
     $     ,r2z,r3x,r3y,r3z, vp1,vp2,vp3,v31,dv1,d13,d23,dv2 , theta1
     $     ,theta2,h,vvp,rcb,vorp,rad

      DOUBLE PRECISION spx(sdim),spy(sdim),spz(sdim), nx0(sdim),ny0(sdim
     $     ),nz0(sdim), nx1(sdim),ny1(sdim),nz1(sdim), rtx(sdim)
     $     ,rty(sdim),rtz(sdim), mrt(sdim)

      DOUBLE PRECISION tic(sdim,sdim),bic(sdim,sdim),nwicm(sdim,sdim),
     $     bcx1(rdim,sdim),bcy1(rdim,sdim),tw(rdim,sdim), twt,t0r(rdim)

      double precision bcx(rdim,sdim), bcy(rdim,sdim), bcz(rdim,sdim)
      DOUBLE PRECISION bicm(sdim,sdim)

      DOUBLE PRECISION dp

      INTEGER lengthnw(sdim)

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst,ycomp, zcomp, pitch(sdim),pit

      external vorp,planf


      common/pidata/pi
      save/pidata/

      common/nearwake/nx1,ny1,nz1
      save/nearwake/

      dpn=20.0D0*pi/180.0D0
      rcn=rcb*crd

      do r=1,nr
         if (abs(twt) .gt. 100) then
            do s=1,ns+1
               pitch(s)=t0r(r)*(rad/spx(s)-1.0D+0)
            enddo
         else
            do s=1,ns+1
               pitch(s)=tw(r,s)
            enddo
         endif
         do s=1,ns
            pit=0.5D+0*(pitch(s)+pitch(s+1))
            bcx(r,s)=bcx1(r,s)
            bcy(r,s)=bcy1(r,s)*cos(pit)
            bcz(r,s)=bcy1(r,s)*sin(pit)
         enddo
      enddo


      do s=1,ns+1

         nx0(s)=spx(s)
         ny0(s)=spy(s)
         nz0(s)=spz(s)

         nx1(s)=spx(s)
         ny1(s)=-spx(s)*tan(dpn)-crd
         ny1(s)=ny1(s)*cos(pitch(s))
         nz1(s)=ny1(s)*sin(pitch(s))
         lengthnw(s)=nint(dpn/dp)
      end do

C-$   Near wake influence coefficients
      r=1
      do j=1,ns+1
         rtx(j)=nx1(j)-nx0(j)
         rty(j)=ny1(j)-ny0(j)
         rtz(j)=nz1(j)-nz0(j)
         mrt(j)=sqrt( rtx(j)**2 + rty(j)**2 + rtz(j)**2 )
      end do

      do i=1,ns
         do j=1,ns+1

            r1x=bcx(r,i)-nx0(j)
            r1y=bcy(r,i)-ny0(j)
            r1z=bcz(r,i)-nz0(j)
            mr1=sqrt( r1x**2 + r1y**2 + r1z**2 )

            r2x=bcx(r,i)-nx1(j)
            r2y=bcy(r,i)-ny1(j)
            r2z=bcz(r,i)-nz1(j)
            mr2=sqrt( r2x**2 + r2y**2 + r2z**2 )

            vp1=-r1y*rtz(j)+r1z*rty(j)
            vp2=+r1x*rtz(j)-r1z*rtx(j)
            vp3=-r1x*rty(j)+r1y*rtx(j)
            v31=sqrt( vp1**2 + vp2**2 + vp3**2 )

            d13=r1x*rtx(j)+r1y*rty(j)+r1z*rtz(j)
            d23=r2x*rtx(j)+r2y*rty(j)+r2z*rtz(j)

            theta1=acos(d13/(mr1*mrt(j)))
            theta2=acos(d23/(mr2*mrt(j)))

            h=mr1*sin(theta1)

            vvp=vorp(h,rcn*planf(bcx(r,i)/rad,taperst,taper))
            tic(i,j)=vvp*(cos(theta1)-cos(theta2))
            ycomp=tic(i,j)*vp2/v31
            zcomp=tic(i,j)*vp3/v31
            tic(i,j)=zcomp*cos(pitch(i))-ycomp*sin(pitch(i))

         end do
      end do

c >>  Bound circulation influence coefficients: <<
      do r=1,nr
         do i=1,ns
            do j=1,ns

               r1x=bcx(r,i)-spx(j)
               r1y=bcy(r,i)-spy(j)
               r1z=bcz(r,i)-spz(j)
               mr1=sqrt( r1x**2 + r1y**2 + r1z**2 )

               r2x=bcx(r,i)-spx(j+1)
               r2y=bcy(r,i)-spy(j+1)
               r2z=bcz(r,i)-spz(j+1)
               mr2=sqrt( r2x**2 + r2y**2 + r2z**2 )

               r3x=spx(j)-spx(j+1)
               r3y=spy(j)-spy(j+1)
               r3z=spz(j)-spz(j+1)
               mr3=sqrt( r3x**2 + r3y**2 + r3z**2 )

               vp1=-r1y*r3z+r1z*r3y
               vp2=+r1x*r3z-r1z*r3x
               vp3=-r1x*r3y+r1y*r3x
               v31=sqrt( vp1**2 + vp2**2 + vp3**2 )

               dv1=r3x*r1x+r3y*r1y+r3z*r1z
               dv2=r3x*r2x+r3y*r2y+r3z*r2z

               theta1=acos(dv1/(mr3*mr1))
               theta2=acos(dv2/(mr3*mr2))

               h=mr1*sin(theta1)

               if (h .le. 1.D-6) then
                  bic(i,j)=0.D+0
               else
                  vvp=vorp(h,0.D+0*rcb)
                  bic(i,j)=vvp*(cos(theta1)-cos(theta2))
                  ycomp=bic(i,j)*vp2/v31
                  zcomp=bic(i,j)*vp3/v31
                  bic(i,j)=zcomp*cos(pitch(i))-ycomp*sin(pitch(i))
               endif

            end do
         end do
      end do

c     Near-wake influence coefficient matrix -- bound+trailed:
      do i=1,ns
         do j=1,ns
            bicm(i,j)=- bic(i,j) 
            nwicm(i,j)=-( (tic(i,j+1)-tic(i,j)) )
         end do
      end do

      return
      end
