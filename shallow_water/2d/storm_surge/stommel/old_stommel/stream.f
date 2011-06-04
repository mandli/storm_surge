
      double precision function stream(x,y)
      implicit real*8(a-h,o-z)
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /stommel/ dlam,bb,dd,capf,rr,ff,gamma,alfs,
     &		       capa,capb,pp,qq,dgbp2

c     # stream function for Stommel Gyre
c
      pi = datan(1.d0)*4.d0
      u0 = pi
      omega = 0.5d0
      width = 0.5d0

      xbar = alf*(x-x0) + beta*(y-y0)
      ybar = -beta*(x-x0) + alf*(y-y0)
      if (xbar .gt. 0.d0) then
	  stream =  -dgbp2 * dsin(pi*ybar/bb) * 
     &       (pp*dexp(capa*xbar) + qq*dexp(capb*xbar) - 1.d0)
       else
	  stream = 0.d0
       endif
      
      return
      end

