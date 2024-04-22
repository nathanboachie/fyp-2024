C-$	$Id: trimtest.f,v 1.11.2.3 2005/06/11 15:57:37 shreyas Exp $	

      subroutine trimtest(nr,np,ns,nb,nctrl,ip,dp,den,flph,om,rad,
     $     compress, nubeta,bmass,ibeta,betap,kbeta,crd,clairfoil,claoa
     $     ,cttol,fltol, seg,tw,bicm,nwicm,bcx,bcy,twt,cdairfoil,kinv,
     $     sonic,aoai,aoain,jac,Cnni,Cnqi,vr,vry,vrz,gb,anb,cpx, taperst
     $     ,taper,mu,cpz,asr,rr_x,rr_z, t0r,t1c,t1s,b0r,b1c,b1s,ct,cqi
     $     ,cqp,WantTrim,pif,zif,Itern,alpha0)

C-$   Wind tunnel trim. 

      IMPLICIT none

      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'jdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'


      INTEGER r,nr,p,np,s,ns,nb,j,c,nctrl

      DOUBLE PRECISION aza,ip(rdim),dp,aoac,den,flph,om,rad,nubeta,
     $     sonic,kinv,mu, bmass,ibeta,betap,kbeta,crd,cttol,fltol,pi,rtd
     $     ,anb

      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim), b0r(rdim),b1c(rdim
     $     ),b1s(rdim), b0r_0(rdim),b1c_0(rdim),b1s_0(rdim), ct(rdim)
     $     ,ct_0(rdim),ect(rdim),efl(rdim), err(jdim),crr(jdim),seg(sdim
     $     ),ct0,cp0

      DOUBLE PRECISION er0(rdim,jdim),cr0(rdim,jdim), tw(rdim,sdim)
     $     ,nwicm(sdim,sdim),bcx(rdim,sdim), jac(jdim,jdim),bcy(rdim
     $     ,sdim)

      DOUBLE PRECISION bicm(sdim,sdim)

      DOUBLE PRECISION aoae(rdim,adim,sdim),aoai(rdim,adim,sdim),
     $     wlfv(rdim,adim,sdim),vr(rdim,adim,sdim), Cnni(rdim,adim,sdim)
     $     ,Cnqi(rdim,adim,sdim), gb(rdim,adim,sdim),cpx(rdim,adim,sdim)

      CHARACTER*1 trimmed,WantTrim

      DOUBLE PRECISION twt

      INTEGER pif,zif

      DOUBLE PRECISION vry(rdim,adim,sdim), vrz(rdim,adim,sdim)
      DOUBLE PRECISION cqi(rdim)
      DOUBLE PRECISION cqp(rdim)
      DOUBLE PRECISION aflp,apit
      DOUBLE PRECISION aoain(rdim,adim,sdim),alpha0

      INTEGER s75

      DOUBLE PRECISION clairfoil(rdim,adim,sdim), cdairfoil(rdim,adim
     $     ,sdim), claoa(rdim,adim,sdim),compress(rdim,adim,sdim)
      DOUBLE PRECISION Mach(rdim,adim,sdim),Reynolds(rdim,adim,sdim)

      DOUBLE PRECISION planf
      DOUBLE PRECISION taper,taperst

      DOUBLE PRECISION Timex,maxChange

      DOUBLE PRECISION gbold(rdim,adim,sdim)
      DOUBLE PRECISION t0rold(rdim),t1cold(rdim),t1sold(rdim)

      DOUBLE PRECISION summt0,sumt1c,sumt1s,summMax,cht0,cht1c,cht1s
      INTEGER n,itmp
      LOGICAL DontTrim

      INTEGER Itern,count

      DOUBLE PRECISION cpz(rdim,adim,sdim),asr(rdim)
      DOUBLE PRECISION rr_x(rdim),rr_z(rdim)

      DOUBLE PRECISION fmind,fmt


      double precision ahdot(rdim,adim,sdim),alpdot(rdim,adim,sdim)
      double precision ualpha(rdim,adim,sdim),delalp(rdim,adim,sdim)
      double precision xcoeff(rdim,adim,sdim),ycoeff(rdim,adim,sdim)
      double precision uacirc(rdim,adim,sdim),uancir(rdim,adim,sdim)
      double precision secaoa(rdim,adim,sdim)

      external planf


      common/pidata/pi
      save/pidata/

c$$$ get the required thrust thru the common block

      common/ctdata/ct0
      save/ctdata/

      common/unsdata/ualpha,delalp,ahdot,alpdot,xcoeff,ycoeff,uacirc,
     $     secaoa
      save/unsdata/

      cp0=0.00108825D+0
      do r=1,nr
         t0rold(r)=t0r(r)
         t1cold(r)=t1c(r)
         t1sold(r)=t1s(r)
         
         b0r_0(r)=b0r(r)
         b1c_0(r)=b1c(r)
         b1s_0(r)=b1s(r)
         
         do p=1,np
            do s=1,ns
               gbold(r,p,s)=gb(r,p,s)
            end do
         end do
      end do

      rtd=180.D+0/pi
      maxChange=0.D+0
      DontTrim=.false.



      if(WantTrim.eq.'y') 
     &write(*,*)'----------------------- TRIMMING ----------------------
     &--'
      if(WantTrim.eq.'r') 
     &write(*,*)'-------------------- FLAP RESPONSE --------------------
     &--'
      if(WantTrim.eq.'n') 
     &write(*,*)'-------------------- TRIM SETTING  --------------------
     &--'
      write(*,*)'R   C_T      C_Qi      C_Qp     T_0     T_1c    T_1s'
     &,'    B_0     B_1c    B_1s'


      do r=1,nr
         write(*,'(i2,f8.4,1x,f9.7,1x,f9.7,6f8.4)')r,ct(r),cqi(r),cqp(r)
     $        ,t0r(r)*rtd,t1c(r)*rtd,t1s(r)*rtd,b0r(r)*rtd,b1c(r)*rtd
     $        ,b1s(r)*rtd
      end do
      
      if(WantTrim.ne.'n') then  !dont do anything if trm 'n'
         count=0
         trimmed='n'
         
         do while(trimmed .eq. 'n' .and. count .le. 40)
            do r=1,nr
               do p=1,np
                  aza=DBLE(p-1)*dp/DBLE(pif)
                  aoac=t0r(r)+t1c(r)*cos(aza)+t1s(r)*sin(aza)
                  do s=1,ns
                     
                     if(abs(twt).gt.100) aoac=t0r(r)/(bcx(r,s)/rad)
     $                    +t1c(r)*cos(aza)+t1s(r)*sin(aza)
                  
                     aflp=-bcx(r,s)*om*( -b1c(r)*sin(aza)+b1s(r)
     $                    *cos(aza) )/(-vry(r,p,s))
                     aflp=atan(aflp)
                     apit=-om*(-t1cold(r)*sin(aza)+t1sold(r)*cos(aza))
                     apit=atan(apit*bcy(r,s)/(-vry(r,p,s)))
                     aoae(r,p,s)=aoac+tw(r,s) !+aflp+apit
                     aoae(r,p,s)=aoae(r,p,s)+aoain(r,p,s)
                     aoae(r,p,s)=aoae(r,p,s)-xcoeff(r,p,s)-ycoeff(r,p,s)
                  
                     Mach(r,p,s)=abs(vry(r,p,s))/sonic
                     Reynolds(r,p,s)=den*vr(r,p,s)
     &                    *crd*planf(bcx(r,s)/rad,taperst,taper)/kinv
                  end do        !s
               end do           !p
            end do              !r
            
            do r=1,nr
               call airfoil(bcx,aoae,Mach,Reynolds,nr,np,ns,clairfoil,
     $              cdairfoil,claoa,compress,0,alpha0)
            
               call rsp(r,np,ns,nb,ip,dp,den,flph,om,rad,crd,clairfoil,
     $              claoa, nubeta,bmass,ibeta,betap,kbeta,seg,bcx,vr,vry
     $              , vrz,gb,Cnni, Cnqi,sonic,cdairfoil,aoae,aoai,aoain,
     $              compress,bicm,nwicm, mu,Mach,Reynolds,taperst,taper,
     $              b0r_0,b1c_0,b1s_0,ct_0,cqi, cqp,pif,zif,Itern)
            
            end do              !r


            do r=1,nr
               if (nr .eq. 1) then
                  er0(r,1)=ct0-ct_0(r)
               else if (r .eq. 1) then 
                  er0(r,1)=ct0-(ct_0(nr-1)+ct_0(nr))
               else
                  er0(r,1)=(cqi(nr-1)+cqp(nr-1))-(cqi(nr)+cqp(nr))
               endif
c$$$               er0(r,1)=ct0-ct_0(r)
               er0(r,2)=b1c(r)-b1c_0(r)
               er0(r,3)=b1c(r)-b1s_0(r)
               ect(r)=abs( er0(r,1)/ct0 )
               efl(r)=sqrt( er0(r,2)**2 + er0(r,3)**2 )
            end do


            do r=1,nr
c$$$               if((WantTrim .eq. 'y') .and. ((ect(1).gt.cttol).or.
c$$$     $              (efl(r).gt.fltol))) then
               if((WantTrim .eq. 'y') .and. (ect(1).gt.cttol)) then
                  trimmed='n'
               else
                  trimmed='y'
                  write(*,*) r, er0(r,1)
               endif

               write(*,'(i2,f8.4,1x,f9.7,1x,f9.7,6f8.4)')r,ct_0(r),cqi(r
     $              ),cqp(r),t0r(r)*rtd,t1c(r)*rtd,t1s(r)*rtd, b0r_0(r)
     $              *rtd,b1c_0(r)*rtd,b1s_0(r)*rtd
            end do
            
            if (trimmed .eq. 'n') then
               count=count+1
               j=0
               do r=1,nr
                  do c=1,nctrl
                     j=j+1
                     err(j)=er0(r,c)
                  end do
               end do
               
               do r=1,nr
                  call ctrlsolve(r,nr,nctrl,jac,err,crr)
               end do
               
               j=0
               do r=1,nr
                  do c=1,nctrl
                     j=j+1
                     cr0(r,c)=crr(j)
                  end do
               end do
               
               do r=1,nr
                  t0r(r)=t0r(r)+cr0(r,1)
                  t1c(r)=t1c(r)+cr0(r,2)
                  t1s(r)=t1s(r)+cr0(r,3)
               end do
            endif
         enddo
         if(trimmed .eq. 'n') write(*,*) "====> Trim not complete"

         if( (WantTrim.eq.'y') .OR. (WantTrim.eq.'r')  ) then
            do r=1,nr
               b0r(r)=b0r_0(r)
               if (abs(b1c_0(r)) .gt. 1.0D-4) b1c(r)=b1c_0(r)
               if (abs(b1s_0(r)) .gt. 1.0D-4) b1s(r)=b1s_0(r)
            end do
         endif
      endif
      
      return
      end
      
