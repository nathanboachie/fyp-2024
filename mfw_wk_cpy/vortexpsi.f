C-$	$Id: vortexpsi.f,v 1.1.1.1.6.2 2005/06/11 15:57:37 shreyas Exp $	

      subroutine vortexpsi(nr,np,nw,ns,nzt,nk,dp,anb,rad,cpx,spx,gb,rv,
     &                 gv,p)

C-$   Trailed vortex strenghts... 


      IMPLICIT none

      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER p0,z,nzt,pm

      INTEGER nr,np,nw,ns,nk,r,p,w,s,i,j,k,dnp

      DOUBLE PRECISION rad,dp,anb

      DOUBLE PRECISION spx(sdim)

      DOUBLE PRECISION cpx(rdim,adim,sdim), gb(rdim,adim,sdim),rv(rdim
     $     ,adim,wdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim)

      DOUBLE PRECISION gpeak, gsum1, gsum2
      INTEGER npeak,sw,sw2,msw,n1, nstart,nend


      k=nk
      
      do r=1,nr
         gv(r,p,1,1,k)=0.D+0
         rv(r,p,1)=1.D+0*rad
         npeak=ns
         do s=ns-1,2,-1
            if (gb(r,p,s).ge.gb(r,p,s+1) .AND. gb(r,p,s).gt.gb(r,p,s-1))
     $           then
               npeak=s
               goto 1014
            endif
            if (gb(r,p,s).lt.gb(r,p,s+1) .AND. gb(r,p,s).le.gb(r,p,s-1))
     $           then
               goto 1014
            endif
            
         end do                 !s
 1014  gpeak=gb(r,p,npeak)

       gv(r,p,1,1,k)=1.00D+0*gpeak
       rv(r,p,1)=1.D+0*rad
       if (nw .gt. 1) then 
          do w=1,nw,1
                                !n1= ns+1 - (w-1)
                                !n1=3-w
             nend=ns+1-4*(w-1)
             nstart=ns+1-4*w
             if (nstart .lt. 1) nstart=1
             gsum1=gsum2
             gsum2=0.0D0
             rv(r,p,w)=spx(nend)
c$$$          if(w.ne.1.AND.w.ne.nw) gv(r,p,w,1,k)=-gb(r,p,n1)+gb(r,p,n1-1)
c$$$          if(w.eq.nw) gv(r,p,w,1,k)=-gb(r,p,n1)
c$$$          if(n1.gt.npeak.AND.w.ne.1) gv(r,p,w,1,k)=0.D+0
c          if (w .eq. nw) then
c             gv(r,p,w,1,k)=gpeak
c          endif
             do n1=nstart,(nend-1)
                gsum2=gsum2+gb(r,p,n1)
             enddo
             if( w .gt. 1) then
                gv(r,p,w,1,k)=gb(r,p,nstart)-gb(r,p,nend)
             else
                gv(r,p,w,1,k)=gb(r,p,nstart)
             endif
          end do                !w
       endif
      end do                    !r


      return
      end
