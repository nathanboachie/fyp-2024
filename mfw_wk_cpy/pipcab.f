C-$	$Id: pipcab.f,v 1.1.1.1.6.2 2005/03/23 21:24:38 shreyas Exp $	
c----------------------------------------------------------------------|
c    <Copyright (c) 1997 by Ashish Bagai>
c    <Copyright (c) 2000 by Mahendra Bhagwat> 
c----------------------------------------------------------------------|
c
      subroutine pipcab(nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,ip,dp,dz,
     &                sonic,vinfx,vinfy,vinfz,om,asr,b0r,anb,rad,rcb,gb,
     &                t0r,t1c,t1s,b1c,b1s,
     &                lin,gv,rc0,rv,kinv,dcy,gav,sco,sce,gsc,csc,
     &                rpar,rrn,zrn,rbc,zbc,
     &                avgzbc,lami,rr_x,rr_y,rr_z,spx,spy,spz,
     &                cpx,cpy,cpz,pcx,pcy,pcz,
     &                pox,poy,poz,ppx,ppy,ppz,bvx,bvy,bvz,
     &                eqtim)

C-$   Bagai's Pseudo-implicit predictor corrector method. 

      IMPLICIT none
c
        include 'adimpar.inc'
        include 'gdimpar.inc'
        include 'jdimpar.inc'
        include 'kdimpar.inc'
        include 'rdimpar.inc'
        include 'sdimpar.inc'
        include 'wdimpar.inc'
        include 'zdimpar.inc'
c
      INTEGER nr,np,nw,pw,nz,nzt,nk,ns,nb,pif,zif,
     &        r,p,w,z,n,s
c
      DOUBLE PRECISION ip(rdim),dp,dz,vinfx,vinfy,vinfz,om,anb,rad,
     &       rcb,lin,kinv,dcy,gav,rpar,rrn,zrn,sce,
     &       ze,wage,tt11,tt21,tt31,b0,sonic,pbar,qbar,pi
c
      DOUBLE PRECISION rbc(rdim,wdim),zbc(rdim,wdim),
     &       avgzbc(rdim,wdim)
c
      DOUBLE PRECISION asr(rdim),b0r(rdim),lami(rdim),
     &       rr_x(rdim),rr_y(rdim),rr_z(rdim),
     &       spx(sdim),spy(sdim),spz(sdim)
c
      DOUBLE PRECISION b1c(rdim),b1s(rdim)
c
      DOUBLE PRECISION rv(rdim,adim,wdim),
     &       gb(rdim,adim,sdim),
     &       cpx(rdim,adim,sdim),
     &       cpy(rdim,adim,sdim),
     &       cpz(rdim,adim,sdim),
     &       vibx(rdim,adim,sdim),
     &       viby(rdim,adim,sdim),
     &       vibz(rdim,adim,sdim),
     &       bvx(rdim,adim,sdim),
     &       bvy(rdim,adim,sdim),
     &       bvz(rdim,adim,sdim)
c
      DOUBLE PRECISION gv(rdim,adim,wdim,zdim,kdim),
     &       rc0(rdim,adim,wdim,kdim),
     &       ppx(rdim,adim,wdim,zdim),
     &       ppy(rdim,adim,wdim,zdim),
     &       ppz(rdim,adim,wdim,zdim),
     &       pcx(rdim,adim,wdim,zdim),
     &       pcy(rdim,adim,wdim,zdim),
     &       pcz(rdim,adim,wdim,zdim),
     &       pox(rdim,adim,wdim,zdim),
     &       poy(rdim,adim,wdim,zdim),
     &       poz(rdim,adim,wdim,zdim),
     &       vipx(rdim,adim,wdim,zdim),
     &       vipy(rdim,adim,wdim,zdim),
     &       vipz(rdim,adim,wdim,zdim),
     &       vicx(rdim,adim,wdim,zdim),
     &       vicy(rdim,adim,wdim,zdim),
     &       vicz(rdim,adim,wdim,zdim),
     &       eqtim(rdim,adim,wdim,zdim,kdim)
c
      CHARACTER*1 sco
c
      DOUBLE PRECISION gsc(rdim,adim,wdim,zdim,kdim),
     &       csc(rdim,adim,wdim,zdim,kdim)
c
      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim)
c
      DOUBLE PRECISION rcore(rdim,adim,wdim,zdim)
     &               ,squire(rdim,adim,wdim,zdim)
     &              ,osquire(rdim,adim,wdim,zdim)
      DOUBLE PRECISION trb,zet0
c
      INTEGER Itern
c----------------------------------------------------------------------|
c
      common/pidata/pi
      save/pidata/
c
      common/ssmdata/pbar,qbar
      save/ssmdata/
c
      common/coredata/rcore,squire,osquire
      save/coredata/
c----------------------------------------------------------------------|
c
           do r=1,nr
           do p=1,np
           do w=1,nw
           do z=1,nz
              vipx(r,p,w,z)=0.D+0
              vipy(r,p,w,z)=0.D+0
              vipz(r,p,w,z)=0.D+0
              vicx(r,p,w,z)=0.D+0
              vicy(r,p,w,z)=0.D+0
              vicz(r,p,w,z)=0.D+0
           end do
           end do
           end do
           end do
c
c-mb induced velocities
           call indvab(nr,np,nw,pw,nz,nzt,nk,ns,pif,zif,dp,dz,rad,
     &               om,anb,lin,dcy,gav,kinv,gv,rc0,sco,sce,gsc,csc,
     &               rcb,gb,ip,b0r,asr,rr_x,rr_y,rr_z,
     &                t0r,t1c,t1s,b1c,b1s,
     &               spx,spy,spz,pox,poy,poz,bvx,bvy,bvz,
     &               vipx,vipy,vipz,cpx,cpy,cpz,
     &               eqtim)
c
c-mb predictor
           call predab(nr,np,nw,pw,nz,nzt,nk,dp,dz,pif,zif,rad,om,rbc,
     &            zbc,avgzbc,lami,gv,rr_x,rr_y,rr_z,vinfx,vinfy,vinfz,
     &            pox,poy,poz,vipx,vipy,vipz,ppx,ppy,ppz)

c-mb induced velocities
           call indvab(nr,np,nw,pw,nz,nzt,nk,ns,pif,zif,dp,dz,rad,
     &               om,anb,lin,dcy,gav,kinv,gv,rc0,sco,sce,gsc,csc,
     &               rcb,gb,ip,b0r,asr,rr_x,rr_y,rr_z,
     &                t0r,t1c,t1s,b1c,b1s,
     &               spx,spy,spz,ppx,ppy,ppz,bvx,bvy,bvz,
     &               vicx,vicy,vicz,cpx,cpy,cpz,
     &               eqtim)
c
c-mb corrector
           call corrab(nr,np,nw,pw,nz,nzt,nk,dp,dz,pif,zif,rad,om,rbc,
     &            zbc,avgzbc,lami,gv,rr_x,rr_y,rr_z,vinfx,vinfy,vinfz,
     &            rpar,vipx,vipy,vipz,ppx,ppy,ppz,vicx,vicy,
     &            vicz,pcx,pcy,pcz)

c----------------------------------------------------------------------|
c
      return
      end
