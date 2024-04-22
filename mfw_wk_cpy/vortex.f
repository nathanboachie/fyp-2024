C-$	$Id: vortex.f,v 1.2.2.2 2005/06/11 15:57:37 shreyas Exp $	

      subroutine vortex(nr,np,nw,ns,nzt,nk,dp,anb,rad,cpx,spx,gb,rv, gv
     $     ,periodic)


      IMPLICIT none

      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER p0,z,nzt,periodic

      INTEGER nr,np,nw,ns,nk,r,p,w,s,k

      DOUBLE PRECISION rad,dp,anb

      DOUBLE PRECISION spx(sdim)

      DOUBLE PRECISION cpx(rdim,adim,sdim),
     &       gb(rdim,adim,sdim),rv(rdim,adim,wdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim)

      DOUBLE PRECISION gpeak,gsum1,gsum2
      INTEGER npeak,n1,nstart,nend

      k=nk
      do r=1,nr
         do p=1,np
            gv(r,p,1,1,k)=0.D+0
            rv(r,p,1)=1.D+0*rad
            npeak=ns
            do s=ns-1,2,-1

               if (gb(r,p,s).ge.gb(r,p,s+1) .AND. gb(r,p,s).gt.gb(r,p,s
     $              -1)) then
                  npeak=s
                  goto 1014
               endif
               if (gb(r,p,s).lt.gb(r,p,s+1) .AND. gb(r,p,s).le.gb(r,p,s
     $              -1)) then
                  goto 1014
               endif
       
            end do              !ns-1 ... 2
 1014       gpeak=gb(r,p,npeak)

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
c$$$               if(w.ne.1.AND.w.ne.nw) gv(r,p,w,1,k)=-gb(r,p,n1)+gb(r,p
c$$$     $              ,n1-1)
c$$$               if(w.eq.nw) gv(r,p,w,1,k)=-gb(r,p,n1)
c$$$
c$$$               
c$$$               if(n1.gt.npeak.AND.w.ne.1) gv(r,p,w,1,k)=0.D+0
               do n1=nstart,(nend-1)
                  gsum2=gsum2+gb(r,p,n1)
               enddo
               if( w .gt. 1) then
                  gv(r,p,w,1,k)=gb(r,p,nstart)-gb(r,p,nend)
               else
                  gv(r,p,w,1,k)=gb(r,p,nstart)
               endif
            end do              !w
            endif
         end do                 !p
      end do                    !r
      

      if(periodic.eq.1) then
         do r=1,nr
            do z=2,nzt
               do p=1,np
                  p0=p-1
                  if(p0 .le. 0)p0=p0+np
                  do w=1,nw
                     do k=1,nk
                        gv(r,p,w,z,k)=gv(r,p0,w,z-1,k)
                     end do
                  end do
               end do
            end do
         end do
      endif                     !periodic


      return
      end
