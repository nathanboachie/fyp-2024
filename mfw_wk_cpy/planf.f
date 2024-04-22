C-$	$Id: planf.f,v 1.2.2.2 2005/05/10 16:14:30 shreyas Exp $	

      DOUBLE PRECISION function planf(r,taperst,taper)


      IMPLICIT none

      INTEGER rotgeo

      DOUBLE PRECISION r,taperst,taper

      common/rotdata/rotgeo
      save/rotdata/
      

      if(r.gt.taperst) then
         planf=1.D+0+(taper-1.D+0)/(1.D+0-taperst)*(r -taperst)
      else
         planf=1.0D+0
      endif
         
      return
      end
