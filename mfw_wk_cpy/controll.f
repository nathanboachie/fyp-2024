C-$	$Id: controll.f,v 1.20.2.3 2005/06/11 15:57:37 shreyas Exp $	

      subroutine controll(Iteratn,nb,pblade,dp,pif,om,dnp,np,nr,theta
     $     ,t0r,t1c,t1s,ip,asr,pbar,qbar,tw,vinfx,vinfy,vinfz ,rad)

C-$   Subroutine to read in changes in control inputs and/or fight
C-$   conditions. Useful for simulating maneuvers

      IMPLICIT none

      include 'adimpar.inc'
      include 'rdimpar.inc'
      include 'sdimpar.inc'

      INTEGER count,Iteratn,nb,pblade,t0steps,nblade,Psi,dnp, p,pm1,np
     $     ,r,nr,pif,s

      DOUBLE PRECISION t0r(rdim),t1c(rdim),t1s(rdim)
      DOUBLE PRECISION thincr, vinfx, vinfy, vinfz,rad

      DOUBLE PRECISION timet,pi,thetac,dthetac,dp,om,aza
      DOUBLE PRECISION theta(rdim,adim), tw(rdim,sdim)
      DOUBLE PRECISION asr(rdim),ip(rdim)
      DOUBLE PRECISION pbar,qbar
      DOUBLE PRECISION pbar0,qbar0
      DOUBLE PRECISION freq,timeb, t0trm, t1ctrm, t1strm,offset
      DOUBLE PRECISION offset1, offset2, mass, accold, accnew
      DOUBLE PRECISION ctnew, ctold, thrust, clvold, clvnew
      DOUBLE PRECISION alt, delalt,ct0, ctinp
      DOUBLE PRECISION tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
      DOUBLE PRECISION pold, pnew, qold, qnew, vxb, vyb, vzb,pdot,qdot
      INTEGER flip,flip1


      common/changec/count
      save/changec/

      common/rates/pbar0,qbar0
      save/rates/

      common/pidata/pi
      save/pidata/

C-$   Remember trimmed control inputs for initial='b' runs
      common/controls/t0trm,t1ctrm,t1strm
      save/controls/

      common/inpdata/freq,offset1,offset2,mass, ctinp
      save/inpdata/

C-$   variables which change over iterations
      common/oldvals/ctold,ctnew,accold,clvold,alt,pold,qold,flip,flip1
      save/oldvals/

C-$   Retrieve thrust information from the input files
      common/ctdata/ct0
      save/ctdata/

      common/mandata/pdot,qdot
      save/mandata/

c-mb  store angular rates to be applied later
c-mb  note that time is azimuth Psi and not TIME

      timet=2.D+0*pi*DBLE(Iteratn-1)/DBLE(nb)+dp*DBLE(pblade-1)
      if(timet.eq.0.D+0) then

C-$   remember the intial trimmed flight conditions and other things.
        count=0
        t0trm=t0r(1)
        t1ctrm=t1c(1)
        t1strm=t1s(1)
        pbar0=pbar
        qbar0=qbar
        ctinp=ct0
        ctold=thincr
        ctnew=thincr
        accnew=0.0D+0
        clvnew=0.0D+0
        delalt=0.0D+0
        alt=0.0D+0
c$$$  Read in perturbation values
        open(16,file='freq.input',status='old')
        read(16,*)freq
        read(16,*)offset1
        read(16,*)offset2
        read(16,*)mass
        close(16)

c$$$  Initialize the input variables to proper units
        offset1=offset1*2.0D+0*pi
        offset2=offset2*2.0D+0*pi
        freq=freq*pi/180.D+0
        pold=0.0D+0
        qold=0.0D+0
        flip=0
        flip1=0
      endif
      
      if((timet .gt. offset1) .and. (timet .lt. offset2)) then
         read(58,*) tmp1,tmp2,tmp3,tmp4,tmp5
         do r=1,nr
            t0r(r)=(tmp2+2.50283)*pi/180.0D+0
            t1s(r)=-tmp3*pi/180.0D+0
            t1c(r)=-tmp4*pi/180.0D+0
         enddo
         if (mod(flip,4) .eq. 0) then
            read(59,*) tmp1,tmp2,tmp3,tmp4,tmp5
            vinfx=tmp2*0.3048D+0
            vinfy=-tmp3*0.3048D+0
            vinfz=tmp4*0.3048D+0
            read(70,*) tmp1, tmp2, tmp3, tmp4
            pbar=-tmp2*pi/180.0D+0
            qbar=tmp3*pi/180.0D+0
            pdot=(pbar-pold)*om/(dp*4.0D0)
            qdot=(qbar-qold)*om/(dp*4.0D0)
            pold=pbar
            qold=qbar
         endif
         flip=flip+1   
      endif

      return
      end
