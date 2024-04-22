C-$	$Id: vvind.f,v 1.1.1.1.6.3 2005/06/11 15:57:37 shreyas Exp $	

      subroutine vvind(r0,p0,w0,z0,r1,p1,w1,z1,np,nz,nzt,nk, pif,zif,dz
     $     ,om,anb,gv,rc,rc0,lin, trb,px,py,pz,vix,viy,viz,Nb,factor)

c     Free-vortex/free-vortex induced velocity calculations.


      IMPLICIT none

      include 'adimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER r0,p0,w0,z0,r1,p1,w1,z1,np, nk,pif,zif,chi,chi2,p3,p4,k

      DOUBLE PRECISION dz,om,anb,lin,trb,vix,viy,viz,wage, r1x,r1y,r1z
     $     ,mr1,r2x,r2y,r2z,mr2, lx,ly,lz,ml,ang1,ang2,h,ux,uy,uz, umag
     $     ,ca,cb,cot1,cot2,g1,g2,rc,vvp, ic,ic2,icx,icy,icz,icx2,icy2
     $     ,icz2, vorp

      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim), rc0(rdim,adim,wdim
     $     ,kdim), px(rdim,adim,wdim,zdim), py(rdim,adim,wdim,zdim),
     $     pz(rdim,adim,wdim,zdim)

      DOUBLE PRECISION zet0,deltacg, umaginv
      INTEGER nz,nzt,Nb

      DOUBLE PRECISION factor

      external vorp

      common/velcomp/r1x,r1y,r1z,mr1,r2x,r2y,r2z,mr2
      save/velcomp/

      if (z1 .gt. (1+zif)) then
         r1x=r2x
         r1y=r2y
         r1z=r2z
         mr1=mr2
      else
         r1x=px(r0,p0,w0,z0)-px(r1,p1,w1,z1-zif)
         r1y=py(r0,p0,w0,z0)-py(r1,p1,w1,z1-zif)
         r1z=pz(r0,p0,w0,z0)-pz(r1,p1,w1,z1-zif)
         mr1=sqrt(r1x*r1x+r1y*r1y+r1z*r1z)
      endif

      r2x=px(r0,p0,w0,z0)-px(r1,p1,w1,z1)
      r2y=py(r0,p0,w0,z0)-py(r1,p1,w1,z1)
      r2z=pz(r0,p0,w0,z0)-pz(r1,p1,w1,z1)
      mr2=sqrt(r2x*r2x+r2y*r2y+r2z*r2z)

      lx=px(r1,p1,w1,z1)-px(r1,p1,w1,z1-zif)
      ly=py(r1,p1,w1,z1)-py(r1,p1,w1,z1-zif)
      lz=pz(r1,p1,w1,z1)-pz(r1,p1,w1,z1-zif)
      ml=sqrt(lx*lx+ly*ly+lz*lz)

      ux=ly*r1z-lz*r1y
      uy=lz*r1x-lx*r1z
      uz=lx*r1y-ly*r1x
      umag=sqrt(ux*ux+uy*uy+uz*uz)
      h=umag/ml


      if (umag .gt. 1.0D-8) then
         umaginv=1.0D0/umag
         ang1=(lx*r1x+ly*r1y+lz*r1z)/(ml*mr1)
         ang2=(lx*r2x+ly*r2y+lz*r2z)/(ml*mr2)

         vvp=vorp(h,rc)*gv(r1,p1,w1,(z1-1),1)*(ang1-ang2)*umaginv
         vix=ux*vvp
         viy=uy*vvp
         viz=uz*vvp
      else
         vix=0.0D0
         viy=0.0D0
         viz=0.0D0
      endif

      if(.false.) then
      if ( (mr1 .eq. 0.D+0) .or. (mr2 .eq. 0.D+0) ) then
         ang1=0.D+0
         ang2=0.D+0
         h=0.D+0
      else
         ang1=acos((lx*r1x+ly*r1y+lz*r1z)/(ml*mr1))
         ang2=acos((lx*r2x+ly*r2y+lz*r2z)/(ml*mr2))
         h=mr1*sin(ang1)
      endif

      if(r0.ne.r1 .OR. p0.ne.p1 .OR. w0.ne.w1) then
         if((mr1-3.0D+0*rc).lt.2.0D-15) h=h*00.D+0
         if((mr2-3.0D+0*rc).lt.2.0D-15) h=h*00.D+0
      endif

      ux=ly*r1z-lz*r1y
      uy=lz*r1x-lx*r1z
      uz=lx*r1y-ly*r1x
      umag=sqrt(ux*ux+uy*uy+uz*uz)
      umaginv=1.0D0/umag

      if ( (umag .eq. 0.D+0) .or. (h .eq. 0.D+0) ) then
         ux=0.D+0
         uy=0.D+0
         uz=0.D+0
         ca=0.D+0
         cb=0.D+0
         vix=0.D+0
         viy=0.D+0
         viz=0.D+0
      else
         ux=ux*umaginv
         uy=uy*umaginv
         uz=uz*umaginv

         cot1=1.D+0/tan(ang1)
         cot2=1.D+0/tan(ang2)
         do k=1,nk
            g1=gv(r1,p1,w1,(z1-1),k)
            g2=lin*gv(r1,p1,w1,(z1-0),k) +(1.D+0-lin)*gv(r1,p1,w1,(z1-1)
     $           ,k)
            ca=(g1*cot2-g2*cot1)/(cot2-cot1)
            cb=(g2-g1)/(cot2-cot1)
            vvp=vorp(h,rc)
            ic=vvp*ca*(cos(ang1)-cos(ang2))
            ic2=vvp*cb*(sin(ang2)-sin(ang1))
            icx=ic*ux
            icy=ic*uy
            icz=ic*uz
            icx2=ic2*ux
            icy2=ic2*uy
            icz2=ic2*uz
            vix=icx+icx2
            viy=icy+icy2
            viz=icz+icz2
         end do
      endif
      endif
      return
      end
