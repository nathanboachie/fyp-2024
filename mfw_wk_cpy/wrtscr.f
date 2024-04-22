C-$	$Id: wrtscr.f,v 1.1.1.1.6.1 2005/03/10 18:49:33 shreyas Exp $	

      subroutine wrtscr(nr,nb,nw,pw,ns,nz,np,pif,zif, mu,muc,ft,bct,dp
     $     ,dz,rad,crd,om,ct,asr,lsr,sco, rr_x,rr_y,rr_z)


      IMPLICIT none 

      include 'rdimpar.inc'
      
      INTEGER nr,nb,nw,pw,ns,nz,np,pif,zif,r

      DOUBLE PRECISION mu,muc,ft,bct,dp,dz,rad,crd,om,pi,deg

      DOUBLE PRECISION ct(rdim),asr(rdim), rr_x(rdim),rr_y(rdim)
     $     ,rr_z(rdim), lsr(rdim)

      CHARACTER*1 sco


      common/pidata/pi
      save/pidata/

      deg=180.D+0/pi

      write(*,*)

      write(*,19)' --------------------------------------------------'
      write(*,20)'            MARYLAND FREE-WAKE ANALYSIS            '
      write(*,19)' --------------------------------------------------'

      write(*,*)

      write(*,20)' [         All units SI where applicable          ]'

      write(*,*)

      write(*,21)' Number of Rotors .................... : ',nr
      write(*,21)' Number of Blades .................... : ',nb
      write(*,21)' Number of Blade Segments ............ : ',ns
      write(*,21)' Number of Vortex Trailers ........... : ',nw
      write(*,21)' Number of Free Vortex Trailers ...... : ',nw-pw
      write(*,21)' Number of Collocation Points ........ : ',nz
      write(*,21)' Number of Azimuthal Locations ....... : ',np
      write(*,24)' Shed Circulation Contribution ....... : ',sco
      write(*,21)' Azimuthal Interpolation Factor ...... : ',pif
      write(*,21)' Vortex Filament Interpolation Factor. : ',zif

      write(*,*)

      do r=1,nr

      write(*,23)' Thrust Coefficient ............... R',r,' : ',ct(r)
      write(*,23)' Rotor Long Shaft Angle ........... R',r,' : ',
     &             asr(r)*deg
      write(*,23)' Rotor Lat Shaft Angle ............ R',r,' : ',
     &             lsr(r)*deg
      write(*,23)' Rotor Hub Center (x) ............. R',r,' : ',
     &             rr_x(r)
      write(*,23)' Rotor Hub Center (y) ............. R',r,' : ',
     &             rr_y(r)
      write(*,23)' Rotor Hub Center (z) ............. R',r,' : ', 
     &             rr_z(r)

      write(*,*)

      end do

      write(*,22)' Advance Ratio ....................... : ',mu
      write(*,22)' Climb Speed Ratio ................... : ',muc
      write(*,22)' Number of Free Turns ................ : ',ft
      write(*,22)' Number of B.C. Turns ................ : ',bct
      write(*,22)' Azimuthal Resolution ................ : ',dp*deg
      write(*,22)' Vortex Filament Resolution .......... : ',dz*deg
      write(*,22)' Blade Radius ........................ : ',rad
      write(*,22)' Blade Chord ......................... : ',crd
      write(*,22)' Rotor Rotational Frequency .......... : ',om

      write(*,*)
      write(*,19)' --------------------------------------------------'

19    format(a51)
20    format(a51)
21    format(a41,i10)
22    format(a41,f10.5)
23    format(a37,i1,a3,3(f10.5))
24    format(a41,a10)


      return
      end
