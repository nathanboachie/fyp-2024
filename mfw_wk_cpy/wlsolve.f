C-$	$Id: wlsolve.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine wlsolve(r,p,ns,compress,bicm,nwicm,wlfv, gb)

C-$   Weissinger-L solver

      IMPLICIT none

      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER r,p,s,s1,s2,ns,prow,ip1,ipm(sdim)

      DOUBLE PRECISION maxw,w1,sum,per1

      DOUBLE PRECISION per(sdim),bip(sdim),gb1(sdim)

      DOUBLE PRECISION nwicm(sdim,sdim),w(sdim,sdim)

      DOUBLE PRECISION bicm(sdim,sdim)
      DOUBLE PRECISION compress(rdim,adim,sdim)

      DOUBLE PRECISION wlfv(rdim,adim,sdim), gb(rdim,adim,sdim)

      DOUBLE PRECISION pi

      common/pidata/pi
      save/pidata/

c >>  Pass influence coefficient matrix am(sdim,sdim) into          <<
c >>  working matrix w(sdim,sdim) and initialize permutation vector <<
      do s=1,ns
         ipm(s)=s
         do s1=1,ns
            w(s,s1)=bicm(s,s1)+nwicm(s,s1)*compress(r,p,s)
         end do
      end do

c >>  L-U Factorization  <<
c >>  Permutation vector <<
      do s=1,ns
         per(s)=abs(w(s,1))
         do s1=2,ns
            if (per(s) .lt. abs(w(s,s1))) then
               per(s)=abs(w(s,s1))
            endif
         end do
      end do

c >>  Row with pivot element <<
      do s=1,ns-1
c         prow=0
c         maxw=0.D+0
         prow=1
         maxw=abs(w(1,s))/per(1)         

         do s1=s,ns
            if (maxw .lt. (abs(w(s1,s))/per(s1))) then
               maxw=(abs(w(s1,s)/per(s1)))
               prow=s1
            endif
         end do
c     >>   Interchanging rows <<
         do s1=1,ns
            w1=w(s,s1)
            w(s,s1)=w(prow,s1)
            w(prow,s1)=w1
         end do
c     >>   Interchanging permutation vector <<
         ip1=ipm(s)
         ipm(s)=ipm(prow)
         ipm(prow)=ip1
         per1=per(s)
         per(s)=per(prow)
         per(prow)=per1
c     >>   Replace 0's of U matrix with multipliers <<
         do s1=s+1,ns
            w(s1,s)=w(s1,s)/w(s,s)
            do s2=s+1,ns
               w(s1,s2)=w(s1,s2)-w(s1,s)*w(s,s2)
            end do
         end do
      enddo

c >>  Apply permutation vector to RHS forcing vector      <<
c >>  and solve system of linear equations using forward  <<
c >>  and backward substitution                           <<

      do s=1,ns
         bip(s)=wlfv(r,p,ipm(s))
c     >>    Forward substitution <<
         gb1(s)=bip(s)
         do s1=1,s-1
            gb1(s)=gb1(s)-w(s,s1)*gb1(s1)
         end do
      end do
c     >>   Backward sustitution  <<
      gb(r,p,ns+1)=0.D+0
      do s=ns,1,-1
         sum=0.D+0
         do s1=s+1,ns
            sum=sum+w(s,s1)*gb(r,p,s1)
         end do
         gb(r,p,s)=(1.D+0/w(s,s))*(gb1(s)-sum)
      end do


      return
      end
