C-$	$Id: prej.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine prej(t0r,t1c,t1s,b0r,b1c,b1s,ct,clairfoil,claoa,crd,
     $     sonic,nr,np,ns,nb,ip,dp,den,flph,om,rad,bicm,nwicm,nubeta,
     $     bmass,ibeta,betap,kbeta,seg,bcx,tw,vr,aoai,nctrl,anb,cpx,twt,
     $     taperst,taper,mu,cpz,asr,rr_x,rr_z, jac,pif,zif,cdairfoil,vry
     $     ,vrz,compress,kinv,alpha0)

C-$   Pre-Jacobian -- Initializes the matrix to start off trim routine



      IMPLICIT none


      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'


      INTEGER nr,np,ns,nb,r,r1,p,s,nctrl,c,c1,i,j

      DOUBLE PRECISION ip(rdim),
     &       dp,den,flph,om,rad,nubeta,bmass,ibeta,betap,
     &       kbeta,crd,aza,aoac,dctrl,pi,sonic,anb,kinv,mu

      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim),
     &       b0r(rdim),b1c(rdim),b1s(rdim),
     &       ct(rdim),seg(sdim)

      DOUBLE PRECISION nwicm(sdim,sdim),tw(rdim,sdim),
     &       bcx(rdim,sdim),ctrl(rdim,jdim),rsp0(rdim,jdim),
     &       rsp1(rdim,jdim),jac(jdim,jdim)

      DOUBLE PRECISION bicm(sdim,sdim)

      DOUBLE PRECISION vr(rdim,adim,sdim),
     &       aoai(rdim,adim,sdim),
     &       aoae(rdim,adim,sdim),
     &       wlfv(rdim,adim,sdim),
     &       gb(rdim,adim,sdim),
     &       Cnni(rdim,adim,sdim),Cnqi(rdim,adim,sdim),
     &       cpx(rdim,adim,sdim)

      DOUBLE PRECISION pj(rdim,rdim,3,3)

      DOUBLE PRECISION twt

      INTEGER pif,zif

      DOUBLE PRECISION vry(rdim,adim,sdim),
     &                 vrz(rdim,adim,sdim)
      DOUBLE PRECISION cqi(rdim)
      DOUBLE PRECISION cqp(rdim)

      INTEGER s75

      DOUBLE PRECISION clairfoil(rdim,adim,sdim),
     &       cdairfoil(rdim,adim,sdim),
     &       claoa(rdim,adim,sdim),compress(rdim,adim,sdim),alpha0
      DOUBLE PRECISION Mach(rdim,adim,sdim),Reynolds(rdim,adim,sdim)

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst

      DOUBLE PRECISION cpz(rdim,adim,sdim),asr(rdim)
      DOUBLE PRECISION rr_x(rdim),rr_z(rdim)

      external planf


      common/pidata/pi
      save/pidata/


      dctrl=1.00D+0*pi/180.D+0

      do r=1,nr

         ctrl(r,1)=t0r(r)
         ctrl(r,2)=t1c(r)
         ctrl(r,3)=t1s(r)
         rsp0(r,1)=ct(r)
         rsp0(r,2)=b1c(r)
         rsp0(r,3)=b1s(r)

      end do



      do r=1,nr
         do c=1,nctrl
            ctrl(r,c)=ctrl(r,c)+dctrl

            do r1=1,nr
               do p=1,np
                  aza=DBLE(p-1)*dp/DBLE(pif)
                  aoac=ctrl(r1,1)+ctrl(r1,2)*cos(aza)+ctrl(r1,3)*sin(aza
     $                 )
                  do s=1,ns
                     if(abs(twt).gt.100) aoac=ctrl(r1,1)/(bcx(r,s)/rad)
     $                    +ctrl(r1,2)*cos(aza)+ctrl(r1,3)*sin(aza)

                     aoae(r1,p,s)=aoac+aoai(r1,p,s)+tw(r1,s)

                     Mach(r1,p,s)=abs(vry(r1,p,s))/sonic
                     Reynolds(r1,p,s)=den*vr(r1,p,s)
     &                    *crd*planf(bcx(r1,s)/rad,taperst,taper)/kinv

                  end do        !s
               end do           !p
            end do              !r1


            call airfoil(bcx,aoae,Mach,Reynolds,nr,np,ns, clairfoil
     $           ,cdairfoil,claoa,compress,0,alpha0)


            do r1=1,nr
               call rsp(r1,np,ns,nb,ip,dp,den,flph,om,rad,crd,clairfoil
     $              ,claoa,nubeta,bmass,ibeta,betap,kbeta,seg,bcx,vr,vry
     $              ,vrz,gb,Cnni,Cnqi,sonic,cdairfoil,aoae,aoai,aoai
     $              ,compress,bicm,nwicm,mu,Mach,Reynolds,taperst,taper
     $              ,b0r,b1c,b1s,ct,cqi,cqp,pif,zif,0)

               rsp1(r1,1)=ct(r1)
               rsp1(r1,2)=b1c(r1)
               rsp1(r1,3)=b1s(r1)

               do c1=1,nctrl
                  pj(r,r1,c,c1)=(rsp1(r1,c1)-rsp0(r1,c1))/dctrl
               end do

            end do              !r1


            ctrl(r,c)=ctrl(r,c)-dctrl
         end do
      end do



      i=0


      do r1=1,nr
         do c1=1,nctrl
            i=i+1
            j=0
            do r=1,nr
               do c=1,nctrl
                  j=j+1
                  jac(i,j)=pj(r,r1,c,c1)
               end do
            end do
         end do
      end do

      if (.false.) then
         write(*,*)
         write(*,'(a46)')'  Pre-Trim Jacobian Matrix : '
         write(*,*)


         do r1=1,nr
            do c1=1,nctrl
               write(*,'(3(3E10.3,2x))')((pj(r,r1,c,c1),c=1,nctrl),r=1
     $              ,nr)
            end do
            write(*,*)
         end do
         write(*,*)
      endif

      return
      end
