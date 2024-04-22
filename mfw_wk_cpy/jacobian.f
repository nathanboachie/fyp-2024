C-$	$Id: jacobian.f,v 1.8.2.2 2005/03/23 21:24:38 shreyas Exp $	

      subroutine jacobian(nr,np,nw,pw,nz,ns,nk,nb,nctrl,dp,dz,ip,anb,
     $     rad,crd,clairfoil,claoa,dcy,kinv,om,lin,gvin,rc0,asr,b0r0,
     $     b1c0,b1s0,rr_x,rr_y,rr_z,pcx,pcy,pcz,cpx,cpy,cpz, t0r,t1c,t1s
     $     ,vinfx,vinfy,vinfz,bcx,bcy,lengthnw, tw,bicm,nwicm,den,flph
     $     ,nbe,bmass,ibe,bep,sonic, kbe,seg,spx,spy,spz,sco,sce,gsc,csc
     $     ,vr0,aoai0,aoain0,twt,nzt, taperst,taper, jac,pif,zif,lsr
     $     ,cdairfoil,rcb,gbin,compress,eqtim,alpha0)

C-$   Jacobian matrix for the trim routine

      IMPLICIT none


      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'


      INTEGER nzt

      INTEGER nr,np,nw,pw,nz,ns,nk,nb,r,p,k,s,i,j,c,r1,c1,nctrl

      DOUBLE PRECISION pi,dp,dz,ip(rdim),anb,rad,dcy,kinv,
     &       om,lin,gav,aza,aoac,mu,
     &       crd,dctrl,vinfx,vinfy,vinfz,den,flph,nbe,
     &       bmass,ibe,bep,kbe,sce,sonic

      DOUBLE PRECISION asr(rdim),b0r(rdim),lsr(rdim),
     &       t0r(rdim),t1c(rdim),t1s(rdim),
     &       rr_x(rdim),rr_y(rdim),rr_z(rdim),b1s(rdim),b1c(rdim),
     &       ct(rdim),seg(sdim),spx(sdim),b0r0(rdim),b1c0(rdim),
     &       b1s0(rdim)

      DOUBLE PRECISION spy(sdim),spz(sdim),
     &       bvx(rdim,adim,sdim),
     &       bvy(rdim,adim,sdim),
     &       bvz(rdim,adim,sdim)

      DOUBLE PRECISION bcx(rdim,sdim),nwicm(sdim,sdim),
     &       rsp0(rdim,jdim),
     &       rsp1(rdim,jdim),ctrl(rdim,jdim),jac(jdim,jdim),
     &       tw(rdim,sdim),bcy(rdim,sdim)

      DOUBLE PRECISION bicm(sdim,sdim)

      DOUBLE PRECISION cpx(rdim,adim,sdim),
     &       cpy(rdim,adim,sdim),
     &       cpz(rdim,adim,sdim),
     &       wlfv(rdim,adim,sdim),
     &       aoae(rdim,adim,sdim),
     &       aoanc(rdim,adim,sdim),
     &       aoai(rdim,adim,sdim),
     &       aoai0(rdim,adim,sdim),
     &       aoain(rdim,adim,sdim),
     &       aoain0(rdim,adim,sdim),
     &       aflp(rdim,adim,sdim),
     &       apit(rdim,adim,sdim),
     &       gb(rdim,adim,sdim),
     &       gbin(rdim,adim,sdim),
     &       gb0(rdim,adim,sdim),
     &       rv(rdim,adim,wdim),
     &       vry(rdim,adim,sdim),
     &       vrz(rdim,adim,sdim),
     &       vr(rdim,adim,sdim),
     &       vr0(rdim,adim,sdim),
     &       vibx(rdim,adim,sdim),
     &       viby(rdim,adim,sdim),
     &       vibz(rdim,adim,sdim),
     &       vibx0(rdim,adim,sdim),
     &       viby0(rdim,adim,sdim),
     &       vibz0(rdim,adim,sdim),
     &       vibxi(rdim,adim,sdim),
     &       vibyi(rdim,adim,sdim),
     &       vibzi(rdim,adim,sdim),
     &       vibxi0(rdim,adim,sdim),
     &       vibyi0(rdim,adim,sdim),
     &       vibzi0(rdim,adim,sdim),
     &       Cnni(rdim,adim,sdim),Cnqi(rdim,adim,sdim)

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim),
     &       gvin(rdim,adim,wdim,zdim,kdim),
     &       rc0(rdim,adim,wdim,kdim),
     &       pcx(rdim,adim,wdim,zdim),
     &       pcy(rdim,adim,wdim,zdim),
     &       pcz(rdim,adim,wdim,zdim),
     &       eqtim(rdim,adim,wdim,zdim,kdim)

      DOUBLE PRECISION pj(rdim,rdim,3,3)

      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim),
     &       csc(rdim,adim,wdim,zdim,kdim)

      CHARACTER*1 sco

      DOUBLE PRECISION twt

      INTEGER pif,zif
      DOUBLE PRECISION cqi(rdim)
      DOUBLE PRECISION cqp(rdim)
      DOUBLE PRECISION rcb

      INTEGER s75,w

      DOUBLE PRECISION clairfoil(rdim,adim,sdim),
     &       cdairfoil(rdim,adim,sdim),
     &       claoa(rdim,adim,sdim),compress(rdim,adim,sdim),alpha0
      DOUBLE PRECISION Mach(rdim,adim,sdim),Reynolds(rdim,adim,sdim)

      INTEGER lengthnw(sdim)

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst


      external planf


      common/pidata/pi
      save/pidata/


      mu=vinfx/om/rad
      do r=1,nr
         b0r(r)=b0r0(r)
         b1c(r)=b1c0(r)
         b1s(r)=b1s0(r)
      end do
      
      call rtof(nr,np,ns,pif,dp,ip,b0r0,b1c0,b1s0,asr,lsr 
     $     ,rr_x,rr_y,rr_z,spx,spy,spz, bvx,bvy,bvz)


      do r=1,nr
         gav=0.D+0
         do p=1,np
            do k=1,nk
               gav=gav+abs(gvin(r,p,1,1,k))
            end do
         end do
      end do
      gav=gav/DBLE(nr*np)


      call rindvt(nr,np,nw,pw,nz,ns,nk,dp,dz,ip,anb,crd,rad, dcy,kinv,om
     $     ,lin,gav,gvin,rc0,asr,b0r0,b1c0,b1s0, rr_x,rr_y,rr_z,pcx,pcy
     $     ,pcz,cpx,cpy,cpz, sco,sce,gsc,csc,vibx0,viby0,vibz0, vibxi0
     $     ,vibyi0,vibzi0,compress,lengthnw, pif,zif,bvx,bvy,bvz,rcb
     $     ,gbin,bicm,nwicm,.false., eqtim,lsr)


      do r=1,nr
         do p=1,np
            aza=DBLE(p-1)*dp/DBLE(pif)
            aoac=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)
            do s=1,ns
               if(abs(twt).gt.100) aoac=t0r(r)/(bcx(r,s)/rad)+t1c(r)
     $              *cos(aza)+t1s(r)*sin(aza)
               call aoa(r,p,s,cpx,cpy,cpz,vinfx,vinfy,vinfz,vibx0,viby0
     $              ,vibz0,vibxi0,vibyi0,vibzi0,nr,np,ns,dp,ip,b0r0,b1c0
     $              ,b1s0,t1c,t1s,asr,om,bcx,bcy,aoac,tw,aza,pif,aoai0
     $              ,aoain0,aflp,apit,vr0,vry,vrz,aoae,aoanc,lsr)
               Mach(r,p,s)=abs(vry(r,p,s))/sonic
               Reynolds(r,p,s)=den*vr0(r,p,s) *crd*planf(bcx(r,s)/rad
     $              ,taperst,taper)/kinv
            end do              !s
         end do                 !p
      end do                    !r


      call airfoil(bcx,aoae,Mach,Reynolds,nr,np,ns, clairfoil,cdairfoil
     $     ,claoa,compress,0,alpha0)


      do r=1,nr
         call rsp(r,np,ns,nb,ip,dp,den,flph,om,rad,crd,clairfoil,claoa,
     $        nbe,bmass,ibe,bep,kbe,seg,bcx,vr0,vry,vrz,gb0,Cnni,Cnqi
     $        ,sonic, cdairfoil,aoae,aoai0,aoain0,compress,bicm,nwicm,mu
     $        , Mach,Reynolds, taperst,taper,b0r,b1c,b1s,ct,cqi,cqp,pif
     $        ,zif,0)
      enddo

      do r=1,nr
         ctrl(r,1)=t0r(r)
         ctrl(r,2)=t1c(r)
         ctrl(r,3)=t1s(r)

         if (nr .eq. 1) then
            rsp0(r,1)=ct(r)
         else if (r .eq. 1) then 
            rsp0(r,1)=ct(r)+ct(nr)
         else
            rsp0(r,1)=(cqi(nr)+cqp(nr))-(cqi(nr-1)+cqp(nr-1))
         endif
            
         rsp0(r,2)=b1c(r)
         rsp0(r,3)=b1s(r)

         if (abs(ctrl(r,1)) .le. 1.D-8) ctrl(r,1)=0.D+0
         if (abs(ctrl(r,2)) .le. 1.D-8) ctrl(r,2)=0.D+0
         if (abs(ctrl(r,3)) .le. 1.D-8) ctrl(r,3)=0.D+0
         if (abs(rsp0(r,1)) .le. 1.D-8) rsp0(r,1)=0.D+0
         if (abs(rsp0(r,2)) .le. 1.D-8) rsp0(r,2)=0.D+0
         if (abs(rsp0(r,3)) .le. 1.D-8) rsp0(r,3)=0.D+0
      end do                    !r


      dctrl=0.25D+0*pi/180.D+0

      do r=1,nr
         do c=1,nctrl
            if (abs(ctrl(r,1)) .le. 1.D-8) ctrl(r,1)=0.D+0
            if (abs(ctrl(r,2)) .le. 1.D-8) ctrl(r,2)=0.D+0
            if (abs(ctrl(r,3)) .le. 1.D-8) ctrl(r,3)=0.D+0

            ctrl(r,c)=ctrl(r,c)+dctrl

            do r1=1,nr
               do p=1,np
                  aza=DBLE(p-1)*dp/DBLE(pif) 
                  aoac=ctrl(r1,1)+ctrl(r1,2)*cos(aza)+ctrl(r1,3)*sin(aza
     $                 )
                  do s=1,ns
                     if(abs(twt).gt.100) aoac=ctrl(r1,1)/(bcx(r,s)/rad)
     $                    +ctrl(r1,2)*cos(aza)+ctrl(r1,3)*sin(aza)

                     call aoa(r1,p,s,cpx,cpy,cpz,vinfx,vinfy,vinfz,vibx0
     $                    ,viby0,vibz0,vibxi0,vibyi0,vibzi0,nr,np,ns,dp
     $                    ,ip,b0r0,b1c0,b1s0,t1c,t1s,asr,om,bcx,bcy,aoac
     $                    ,tw,aza,pif,aoai,aoain,aflp,apit,vr,vry,vrz
     $                    ,aoae,aoanc,lsr)

                     Mach(r1,p,s)=abs(vry(r1,p,s))/sonic
                     Reynolds(r1,p,s)=den*vr(r1,p,s) *crd*planf(bcx(r1,s
     $                    )/rad,taperst,taper)/kinv
                  end do        !s
               end do           !p
            end do              !r1


            call airfoil(bcx,aoae,Mach,Reynolds,nr,np,ns, clairfoil
     $           ,cdairfoil,claoa,compress,0,alpha0)


            do r1=1,nr
               call rsp(r1,np,ns,nb,ip,dp,den,flph,om,rad,crd,clairfoil
     $              ,claoa,nbe,bmass,ibe,bep,kbe,seg,bcx,vr,vry,vrz,gb
     $              ,Cnni,Cnqi,sonic, cdairfoil,aoae,aoai,aoain,compress
     $              ,bicm,nwicm,mu,Mach,Reynolds, taperst,taper,b0r,b1c
     $              ,b1s,ct,cqi,cqp,pif,zif,0)

            enddo

            do r1=1,nr

               if (nr .eq. 1) then
                  rsp1(r1,1)=ct(r1)
               else if (r1 .eq. 1) then 
                  rsp1(r1,1)=ct(r1)+ct(nr)
               else
                  rsp1(r1,1)=(cqi(nr)+cqp(nr))-(cqi(nr-1)+cqp(nr-1))
               endif
               rsp1(r1,2)=b1c(r1)
               rsp1(r1,3)=b1s(r1)
               
               if (abs(rsp1(r1,1)) .le. 1.D-8) rsp1(r1,1)=0.D+0
               if (abs(rsp1(r1,2)) .le. 1.D-8) rsp1(r1,2)=0.D+0
               if (abs(rsp1(r1,3)) .le. 1.D-8) rsp1(r1,3)=0.D+0

               do c1=1,nctrl
                  pj(r,r1,c,c1)=(rsp1(r1,c1)-rsp0(r1,c1))/dctrl
               end do

            end do              !r1


            ctrl(r,c)=ctrl(r,c)-dctrl
         end do                 !c-loop 
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

      if(.true.) then
         write(*,*)
         write(*,'(a46)')'  Coupled Jacobian Matrix :  '
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
