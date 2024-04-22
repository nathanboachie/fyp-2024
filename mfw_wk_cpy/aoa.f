C-$	$Id: aoa.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine aoa(r,p,s,cpx,cpy,cpz,vinfx,vinfy,vinfz,vibx,viby,vibz,
     $     vibxin,vibyin,vibzin, nr,np,ns,dp,ip,b0r,b1c,b1s,t1c,t1s,asr
     $     ,om,bcx,bcy, aoac,tw,aza,pif, aoai,aoain,aflp,apit,vr,vry,vrz
     $     ,aoae,aoanc,lsr)

C-$   Calcuate the sectional angle of attack

      IMPLICIT none

      include 'adimpar.inc'
      include 'gdimpar.inc'
      include 'kdimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'
      include 'wdimpar.inc'
      include 'zdimpar.inc'

      INTEGER r,p,s,nr,np,ns

      DOUBLE PRECISION vinfx,vinfy,vinfz,om,aoac,aza,dp,ip(rdim)

      DOUBLE PRECISION b0r(rdim),b1c(rdim),b1s(rdim),asr(rdim), t1c(rdim
     $     ),t1s(rdim),lsr(rdim)

      DOUBLE PRECISION bcx(rdim,sdim),bcy(rdim,sdim), tw(rdim,sdim)

      DOUBLE PRECISION cpx(rdim,adim,sdim), cpy(rdim,adim,sdim),
     $     cpz(rdim,adim,sdim), vibx(rdim,adim,sdim), viby(rdim,adim
     $     ,sdim), vibz(rdim,adim,sdim), vibxin(rdim,adim,sdim),
     $     vibyin(rdim,adim,sdim), vibzin(rdim,adim,sdim), vr(rdim,adim
     $     ,sdim), vfx(rdim,adim,sdim), vfy(rdim,adim,sdim), vfz(rdim
     $     ,adim,sdim), vrx(rdim,adim,sdim), vry(rdim,adim,sdim),
     $     vrz(rdim,adim,sdim), aoai(rdim,adim,sdim), aoain(rdim,adim
     $     ,sdim), aflp(rdim,adim,sdim), apit(rdim,adim,sdim), aoae(rdim
     $     ,adim,sdim),aoanc(rdim,adim,sdim)

      double precision pbar,qbar

      INTEGER pif
      

      common/ssmdata/pbar,qbar
      save/ssmdata/


      vfx(r,p,s)=vinfx-qbar*cpz(r,p,s)
      vfy(r,p,s)=vinfy+pbar*cpz(r,p,s)
      vfz(r,p,s)=vinfz+qbar*cpx(r,p,s)-pbar*cpy(r,p,s)

      call ftor(r,p,s,dp,ip,b0r,asr,lsr,b1c,b1s,vfx,vfy,vfz,vrx,vry,vrz
     $          ,pif)

      vrx(r,p,s)=vrx(r,p,s)+vibx(r,p,s)
      vry(r,p,s)=vry(r,p,s)+viby(r,p,s)
      vrz(r,p,s)=vrz(r,p,s)+vibz(r,p,s)

      vry(r,p,s)=vry(r,p,s)-om*bcx(r,s)

      vr(r,p,s)=sqrt(vrx(r,p,s)**2+vry(r,p,s)**2+vrz(r,p,s)**2 )

      aoai(r,p,s)=atan(vrz(r,p,s)/(-vry(r,p,s)))

      aoain(r,p,s)=atan( (vrz(r,p,s)-vibz(r,p,s)+vibzin(r,p,s)) /(
     $     -vry(r,p,s)+viby(r,p,s)-vibyin(r,p,s)) )
      
      aflp(r,p,s)=-bcx(r,s)*om*( -b1c(r)*sin(aza) + b1s(r)*cos(aza) )/(
     $     -vry(r,p,s))
      aflp(r,p,s)=atan(aflp(r,p,s))

      apit(r,p,s)=-om*(-t1c(r)*sin(aza)+t1s(r)*cos(aza))
      apit(r,p,s)=atan(apit(r,p,s)*bcy(r,s)/(-vry(r,p,s)))

      aoae(r,p,s)=aoac+aoai(r,p,s)+tw(r,s)+ aflp(r,p,s)+apit(r,p,s)

      return
      end
