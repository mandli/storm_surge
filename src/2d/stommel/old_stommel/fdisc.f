c
c
c
c     =================================================
      function fdisc(x,y)
c     =================================================
      implicit double precision (a-h,o-z)
      common/cdisc/ x0,y0,alf,beta,r0,idisc
c
c     # for computing cell averages for initial data that has a
c     # discontinuity along some curve.  fdisc should be negative to the 
c     # left of the curve and positive to the right
c     # idisc specifies the nature of the discontinuity for two
c     # particular cases (a straight line and circle) but this routine
c     # can be modified for any other curve.
c
      go to (10,20,30,40) idisc
c
   10 continue
c     # straight line through (x0,y0) with normal (alf,beta) pointing 
c     # into left state  (for use in cartfr)
c
      fdisc = (x-x0)*alf + (y-y0)*beta
      fdisc = -fdisc
      return
c
   20 continue
c     # outside circle of radius r0:
      fdisc = (x-x0)**2 + (y-y0)**2 - r0**2
      fdisc = -fdisc
c
   30 continue
c     # inside circle of radius r0:
      fdisc = (x-x0)**2 + (y-y0)**2 - r0**2
c
   40 continue
c     # inside sinusoidal channel:
      pi = datan(1.d0)*4.d0
      omega = 0.5d0
      width = 0.5d0
      stream = y - dcos(2.d0*pi*omega*x)
      stream = - dmin1(width, dmax1(0.d0, stream))
      if (stream .eq. 0.d0 .or. stream .eq. -width) then
	  fdisc = 1.d0
	else
	  fdisc = -1.d0
	endif
c
      return
      end
