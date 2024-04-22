c----------------------------------------------------------------------|
c     Copyright violation - MrB
c----------------------------------------------------------------------|
c
      subroutine rtofpsi(nr,np,ns,pif,dp,ip,beta,asr,lsr,
     &                rr_x,rr_y,rr_z,spx,spy,spz,
     &                bvx,bvy,bvz,p)
c     Rotating to fixed frame transformation.
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
        include 'adimpar.inc'
        include 'kdimpar.inc'
        include 'rdimpar.inc'
        include 'sdimpar.inc'
c
      INTEGER nr,np,ns,pif,r,p,s
      INTEGER rotgeo
c
      DOUBLE PRECISION sn,aza,ip(rdim),dp,
     &       tt11,tt12,tt13,
     &       tt21,tt22,tt23,
     &       tt31,tt32,tt33
c
      DOUBLE PRECISION spx(sdim),spy(sdim),spz(sdim),
     &       rr_x(rdim),rr_y(rdim),rr_z(rdim),
     &       asr(rdim),lsr(rdim),beta(rdim,adim),b0
c
      DOUBLE PRECISION bvx(rdim,adim,sdim),
     &       bvy(rdim,adim,sdim),
     &       bvz(rdim,adim,sdim)
     
      common/rotdata/rotgeo
      save/rotdata/

c----------------------------------------------------------------------
c>>mb why use beta_0 here?
c>>mb shouldnt it be the actual flap angle @ each azimuth??
c>>mb  well, Dr A Bugi used beta_0 for reasons forever unknown
c
      do r=1,nr
         sn=1.0D0
         if(rotgeo.eq.1 .or. rotgeo.eq.2) sn=(-1.D+0)**(r-1)

         aza=( ip(r)+DBLE(p-1)*dp/DBLE(pif) ) - ip(r)

         b0=beta(r,p)

         aza=(aza + ip(r))*sn
c-mb sign of rotation for coord transform
c
         tt11= cos(-asr(r))*cos(aza)*cos(b0)+
     &         sin(-asr(r))*cos(lsr(r))*sin(b0)-
     &         sn*sin(-asr(r))*sin(lsr(r))*sin(aza)*cos(b0)
         tt12=-cos(-asr(r))*sin(aza)
     &        -sn*sin(-asr(r))*sin(lsr(r))*cos(aza)
         tt13=-cos(-asr(r))*cos(aza)*sin(b0)
     &           +sn*sin(-asr(r))*sin(lsr(r))*sin(aza)*sin(b0) 
     &           +sin(-asr(r))*cos(lsr(r))*cos(b0)
     
         tt21= sn*cos(lsr(r))*sin(aza)*cos(b0)+sin(lsr(r))*sin(b0)
         tt22= cos(lsr(r))*cos(aza)*sn
         tt23= -sn*cos(lsr(r))*sin(aza)*sin(b0)+sin(lsr(r))*cos(b0)
            
         tt31= -sin(-asr(r))*cos(aza)*cos(b0)
     &           -sn*cos(-asr(r))*sin(lsr(r))*sin(aza)*cos(b0)
     &           +cos(-asr(r))*cos(lsr(r))*sin(b0)
         tt32= sin(-asr(r))*sin(aza)
     &           -sn*cos(asr(r))*sin(lsr(r))*cos(aza)
         tt33= sin(-asr(r))*cos(aza)*sin(b0)
     &           +sn*cos(asr(r))*sin(lsr(r))*sin(aza)*sin(b0) 
     &           +cos(-asr(r))*cos(lsr(r))*cos(b0)


      do s=1,ns+1

         bvx(r,p,s)=rr_x(r)+tt11*spx(s)+
     &                      tt12*spy(s)*sn+
     &                      tt13*spz(s)
         bvy(r,p,s)=rr_y(r)+tt21*spx(s)+
     &                      tt22*spy(s)*sn+
     &                      tt23*spz(s)
         bvz(r,p,s)=rr_z(r)+tt31*spx(s)+
     &                      tt32*spy(s)*sn+
     &                      tt33*spz(s)

      end do !s

      end do !r
c----------------------------------------------------------------------|
c
      return
      end
