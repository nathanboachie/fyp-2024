C-$	$Id: ctrlsolve.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine ctrlsolve(r,nr,nctrl,jac,err,crr)


c-mb basically matrix solver
c-mb     [J]{c}={e}
c-mb -->  {c}=inv[J]{e}


      IMPLICIT none

      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER r,j,k,s2,nr,nctrl,prow,ip1,ipm(sdim)

      DOUBLE PRECISION maxw,w1,sum,per1

      DOUBLE PRECISION per(sdim),bip(sdim),gb1(sdim), err(jdim),crr(jdim
     $     )

      DOUBLE PRECISION w(sdim,sdim),jac(jdim,jdim)


c >>  Pass influence coefficient matrix am(sdim,sdim) into          <<
c >>  working matrix w(sdim,sdim) and initialize permutation vector <<
      do j=1,nr*nctrl
         ipm(j)=j
         do k=1,nr*nctrl
            w(j,k)=jac(j,k)
         end do
      end do

c >>  L-U Factorization <<

c >>  Permutation vector <<
      do j=1,nr*nctrl
         per(j)=abs(w(j,1))
         do k=2,nr*nctrl
            if (per(j) .lt. abs(w(j,k))) then
               per(j)=abs(w(j,k))
            endif
         end do
      end do

c >>  Row with pivot element <<
      do j=1,nr*nctrl-1
         prow=0
         maxw=0.D+0
         do k=j,nr*nctrl
            if (maxw .lt. (abs(w(k,j))/per(k))) then
               maxw=(abs(w(k,j)/per(k)))
               prow=k
            endif
         end do
c >>   Interchanging rows <<
       do k=1,nr*nctrl
          w1=w(j,k)
          w(j,k)=w(prow,k)
          w(prow,k)=w1
       end do
c     >>   Interchanging permutation vector <<
       ip1=ipm(j)
       ipm(j)=ipm(prow)
       ipm(prow)=ip1
       per1=per(j)
       per(j)=per(prow)
       per(prow)=per1
c     >>   Replace 0's of U matrix with multipliers <<
       do k=j+1,nr*nctrl
          w(k,j)=w(k,j)/w(j,j)
          do s2=j+1,nr*nctrl
             w(k,s2)=w(k,s2)-w(k,j)*w(j,s2)
          end do
       end do
      enddo

c >>  Apply permutation vector to RHS forcing vector      <<
c >>  and solve system of linear equations using forward  <<
c >>  and backward substitution                           <<
      do j=1,nr*nctrl
         bip(j)=err(ipm(j))
c     >>    Forward substitution <<
         gb1(j)=bip(j)
         do k=1,j-1
            gb1(j)=gb1(j)-w(j,k)*gb1(k)
         end do
      end do
c     >>   Backward sustitution <<
      crr(nr*nctrl+1)=0.D+0
      do j=nr*nctrl,1,-1
         sum=0.D+0
         do k=j+1,nr*nctrl
            sum=sum+w(j,k)*crr(k)
         end do
         crr(j)=(1.D+0/w(j,j))*(gb1(j)-sum)
      end do


      return
      end
