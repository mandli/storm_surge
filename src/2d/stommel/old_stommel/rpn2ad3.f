
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann-solver for the advection equation
c     #    q_t  +  (u*q)_x + (v*q)_y = 0
c     # where u and v are a given velocity field.
c
c       -----------------------------------------------------------
c     # In conservative form, with interface velocities specified in
c     # the auxiliary variable   
c     # aux(i,j,1)  =  u-velocity at left edge of cell (i,j)
c     # aux(i,j,2)  =  v-velocity at bottom edge of cell (i,j)
c       -----------------------------------------------------------
c
c     # solve Riemann problems along one slice of data.
c     # This data is along a slice in the x-direction if ixy=1
c     #                            or the y-direction if ixy=2.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, wave contains the waves, s the speeds, 
c     # and amdq, apdq the left-going and right-going flux differences,
c     # respectively.  Note that in this advective form, the sum of
c     # amdq and apdq is not equal to a difference of fluxes except in the
c     # case of constant velocities.
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit real*8(a-h,o-z)
c
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension auxl(1-mbc:maxm+mbc, 2)
      dimension auxr(1-mbc:maxm+mbc, 2)
      dimension apdq(1-mbc:maxm+mbc, meqn)
      dimension amdq(1-mbc:maxm+mbc, meqn)
c
c
c     # Set wave, speed, and flux differences:
c     ------------------------------------------
c
c     # set fim1 = f(q_{i-1}) at left boundary for flux differencing below:
      i = 1-mbc
      fim1 = 0.5d0*(qr(i,1) + ql(i,1)) * 
     &            (dmin1(auxl(i+1,ixy),0.d0) 
     &             + dmax1(auxl(i,ixy),0.d0))
c
      do 30 i = 2-mbc, mx+mbc-1
	 u = auxl(i,ixy)
         wave(i,1,1) = ql(i,1) - qr(i-1,1)
	 s(i,1) = u
c
c        # conservative form
c          -----------------
c        # amdq and apdq are chosen as flux differences for the
c        # conservative equation  q_t + (u*q)_x = 0
c
c        # compute the flux at the interface between cells i-1 and i:
	 if (u.gt.0.d0) then
	       f0 = u*qr(i-1,1)
	     else
	       f0 = u*ql(i,1)
	     endif
c
c        # compute a value for the flux in cell i:
c        # note that we have velocities only at the interfaces
	 fi = 0.5d0*(ql(i,1) + qr(i,1)) *
     &             (dmin1(auxl(i+1,ixy),0.d0) + dmax1(u,0.d0))
c
c        # flux differences:
	 amdq(i,1) = f0 - fim1
	 apdq(i,1) = fi - f0
c
	 fim1 = fi
   30    continue
c
      return
      end
