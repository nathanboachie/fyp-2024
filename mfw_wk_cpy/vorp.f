c----------------------------------------------------------------------|
      DOUBLE PRECISION function vorp(x,y)
c----------------------------------------------------------------------|
c Vortex tangential velocity profile model -- 
c n1 (Scully), n2, ninf (Rankine)
c----------------------------------------------------------------------|
c
      IMPLICIT none
c
      DOUBLE PRECISION x,y,pi,vn 
c----------------------------------------------------------------------|
c
      common/pidata/pi
      save/pidata/
c
      common/vndata/vn
      save/vndata/
c----------------------------------------------------------------------|
c
c --- n=1  => Scully vortex tangential velocity profile
c --- n=2  => Bagai/Leishman choice
c --- n=50 => Rankine vortex tangential velocity profile
c
      vorp=x/( 4.D+0*pi*( (x**(2*vn)+y**(2*vn))**(1.D+0/vn) ) )
c----------------------------------------------------------------------|
c
      return
      end
