      subroutine setprob
      implicit double precision (a-h,o-z)
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /stommel/ dlam,bb,dd,capf,rr,ff,gamma,alfs,
     &		       capa,capb,pp,qq,dgbp2

      open(unit=7,file='setprob.data',status='old',form='formatted')
c
c
c     # wall:
      x0 = 1.5d8
      y0 = 0.d0
      read(7,*) x0,y0,theta
      alf = dcos(theta)
      beta = dsin(theta)
c
c     # Stommel parameters:
c
      pi = 4.d0*datan(1.d0)
      pi2 = 2.d0*pi
      dlam = 1.d9
      bb = pi2*1.d8
      dd = 2.d4
      capf = 1.d0
      rr = 0.02d0
      ff = 1.d-13
c     ff = 1.d-13
c     ff = 0.d0
      gamma = capf*pi / (rr*bb)
      alfs = (dd/rr) * ff
      capa = -alfs/2.d0 + dsqrt(alfs**2/4.d0 + (pi/bb)**2)
      capb = -alfs/2.d0 - dsqrt(alfs**2/4.d0 + (pi/bb)**2)
      pp = (1.d0 - dexp(capb*dlam)) / (dexp(capa*dlam) -
     &       dexp(capb*dlam))
      qq = 1.d0 - pp
      dgbp2 = gamma * (bb/pi)**2
      return
      end
