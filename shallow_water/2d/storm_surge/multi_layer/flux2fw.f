c
c
c     =====================================================
      subroutine flux2(ixy,maxm,meqn,maux,mbc,mx,
     &                 q1d,dtdx1d,aux1,aux2,aux3,
     &                 faddm,faddp,gaddm,gaddp,cfl1d,fwave,s,
     &                 amdq,apdq,cqxx,bmasdq,bpasdq,rpn2,rpt2)
c     =====================================================
c
c     # clawpack routine ...  modified for AMRCLAW
c
c--------------------------------------------------------------------
c     # flux2fw is a modified version of flux2 to use fwave instead of wave.
c     # A modified Riemann solver rp2n must be used in conjunction with this
c     # routine, which returns fwave's instead of wave's.
c     # See http://amath.washington.edu/~claw/fwave.html
c
c     # Limiters are applied to the fwave's, and the only significant
c     # modification of this code is in the "do 119" loop, for the
c     # second order corrections.
c
c--------------------------------------------------------------------
c
c
c     # Compute the modification to fluxes f and g that are generated by
c     # all interfaces along a 1D slice of the 2D grid.
c     #    ixy = 1  if it is a slice in x
c     #          2  if it is a slice in y
c     # This value is passed into the Riemann solvers. The flux modifications
c     # go into the arrays fadd and gadd.  The notation is written assuming
c     # we are solving along a 1D slice in the x-direction.
c
c     # fadd(i,.) modifies F to the left of cell i
c     # gadd(i,.,1) modifies G below cell i
c     # gadd(i,.,2) modifies G above cell i
c
c     # The method used is specified by method(2:3):
c
c         method(2) = 1 if only first order increment waves are to be used.
c                   = 2 if second order correction terms are to be added, with
c                       a flux limiter as specified by mthlim.
c
c         method(3) = 0 if no transverse propagation is to be applied.
c                       Increment and perhaps correction waves are propagated
c                       normal to the interface.
c                   = 1 if transverse propagation of increment waves
c                       (but not correction waves, if any) is to be applied.
c                   = 2 if transverse propagation of correction waves is also
c                       to be included.
c
c     Note that if mcapa>0 then the capa array comes into the second
c     order correction terms, and is already included in dtdx1d:
c     If ixy = 1 then
c        dtdx1d(i) = dt/dx                      if mcapa= 0
c                  = dt/(dx*aux(i,jcom,mcapa))  if mcapa = 1
c     If ixy = 2 then
c        dtdx1d(j) = dt/dy                      if mcapa = 0
c                  = dt/(dy*aux(icom,j,mcapa))  if mcapa = 1
c
c     Notation:
c        The jump in q (q1d(i,:)-q1d(i-1,:))  is split by rpn2 into
c            amdq =  the left-going flux difference  A^- Delta q
c            apdq = the right-going flux difference  A^+ Delta q
c        Each of these is split by rpt2 into
c            bmasdq = the down-going transverse flux difference B^- A^* Delta q
c            bpasdq =   the up-going transverse flux difference B^+ A^* Delta q
c        where A^* represents either A^- or A^+.
c

c      #  modifications for GeoClaw

c--------------------------flux2fw_geo.f--------------------------
c     This version of flux2fw.f is modified slightly to be used with
c     step2_geo.f  The only modification is for the first-order
c     mass fluxes, faddm(i,j,1) and faddp(i,j,1), so that those terms are true
c     interface fluxes.
c
c     The only change is in loop 40
c     to revert to the original version, set relimit = .false.
c---------------------last modified 1/04/05-----------------------------

c--------------------------flux2fw.f------------------------------------
c     Changed limiting so that it does not limit if near a dry state in
c     the bottom layer.  This accomplished by not adding in the 2nd 
c     order corrections if near a dry state. (line )
c---------------------last modified 4/28/11-----------------------------

      use multilayer_module, only: rho, dry_limit
      use geoclaw_module

      implicit double precision (a-h,o-z)
      include "call.i"

      external rpn2, rpt2
      dimension  q1d(1-mbc:maxm+mbc, meqn)
      dimension  amdq(1-mbc:maxm+mbc, meqn)
      dimension  apdq(1-mbc:maxm+mbc, meqn)
      dimension  bmasdq(1-mbc:maxm+mbc, meqn)
      dimension  bpasdq(1-mbc:maxm+mbc, meqn)
      dimension  cqxx(1-mbc:maxm+mbc, meqn)
      dimension  faddm(1-mbc:maxm+mbc, meqn)
      dimension  faddp(1-mbc:maxm+mbc, meqn)
      dimension  gaddm(1-mbc:maxm+mbc, meqn, 2)
      dimension  gaddp(1-mbc:maxm+mbc, meqn, 2)
      dimension  dtdx1d(1-mbc:maxm+mbc)
      dimension  aux1(1-mbc:maxm+mbc, maux)
      dimension  aux2(1-mbc:maxm+mbc, maux)
      dimension  aux3(1-mbc:maxm+mbc, maux)
c
      dimension  s(1-mbc:maxm+mbc, mwaves)
      dimension  fwave(1-mbc:maxm+mbc, meqn, mwaves)
c
      logical limit, relimit, dry_l, dry_r
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom



      relimit = .false.
c
      limit = .false.
      do 5 mw=1,mwaves
         if (mthlim(mw) .gt. 0) limit = .true.
   5     continue
c
c     # initialize flux increments:
c     -----------------------------
c
      do 30 jside=1,2
         do 20 m=1,meqn
            do 10 i = 1-mbc, mx+mbc
               faddm(i,m) = 0.d0
               faddp(i,m) = 0.d0
               gaddm(i,m,jside) = 0.d0
               gaddp(i,m,jside) = 0.d0
   10          continue
   20       continue
   30    continue
c
c
c     # solve Riemann problem at each interface and compute Godunov updates
c     ---------------------------------------------------------------------
c
      call rpn2(ixy,maxm,meqn,mwaves,mbc,mx,q1d,q1d,
     &          aux2,aux2,fwave,s,amdq,apdq)
c
c   # Set fadd for the donor-cell upwind method (Godunov)
      if (ixy.eq.1) mu=2
      if (ixy.eq.2) mu=3
      do 40 i=1-mbc+1,mx+mbc-1
         if (icoordsys.eq.2) then
      	  if (ixy.eq.1) dxdc=Rearth*pi/180.d0
	        if (ixy.eq.2) dxdc=Rearth*pi*cos(aux2(i,3))/180.d0
	      else
	       dxdc=1.d0
	      endif

         do m=1,meqn
            faddp(i,m) = faddp(i,m) - apdq(i,m)
            faddm(i,m) = faddm(i,m) + amdq(i,m)
         enddo
         if (relimit) then
            faddp(i,1) = faddp(i,1) + dxdc*q1d(i,mu)
            faddm(i,1) = faddp(i,1)
         endif
   40       continue
c
c     # compute maximum wave speed for checking Courant number:
      cfl1d = 0.d0
      do 50 mw=1,mwaves
         do 50 i=1,mx+1
c          # if s>0 use dtdx1d(i) to compute CFL,
c          # if s<0 use dtdx1d(i-1) to compute CFL:
            cfl1d = dmax1(cfl1d, dtdx1d(i)*s(i,mw),
     &                          -dtdx1d(i-1)*s(i,mw))

   50       continue
c
      if (method(2).eq.1) go to 130
c
c     # modify F fluxes for second order q_{xx} correction terms:
c     -----------------------------------------------------------
c
c     # apply limiter to fwaves:
      if (limit) call limiter(maxm,meqn,mwaves,mbc,mx,fwave,s,mthlim)
c
      do 120 i = 1, mx+1
c
c        # For correction terms below, need average of dtdx in cell
c        # i-1 and i.  Compute these and overwrite dtdx1d:
c
         dtdx1d(i-1) = 0.5d0 * (dtdx1d(i-1) + dtdx1d(i))
c
c        Check to see if we are near a dry state and skip the correction if we
c        are, dry_limit turns this on and off as well
c

         dry_l = q1d(i-1,4) / rho(2) < drytolerance
         dry_r = q1d(i,4) / rho(2) < drytolerance
         if (dry_limit.and.((dry_l.and.(.not.dry_r))
     &           .or.(dry_r.and.(.not.dry_l)))) then
             cycle
         endif
         
         do 120 m=1,meqn
            cqxx(i,m) = 0.d0
            do 119 mw=1,mwaves
c
c              # second order corrections:
               cqxx(i,m) = cqxx(i,m) + dsign(1.d0,s(i,mw))
     &            * (1.d0 - dabs(s(i,mw))*dtdx1d(i-1)) * fwave(i,m,mw)
c
  119          continue
            faddm(i,m) = faddm(i,m) + 0.5d0 * cqxx(i,m)
            faddp(i,m) = faddp(i,m) + 0.5d0 * cqxx(i,m)
  120       continue
c
c
  130  continue
c
       if (method(3).eq.0) go to 999   !# no transverse propagation
c
       if (method(3).eq.2) then
c         # incorporate cqxx into amdq and apdq so that it is split also.
          do 150 i = 1, mx+1
             do 150 m=1,meqn
                amdq(i,m) = amdq(i,m) + cqxx(i,m)
                apdq(i,m) = apdq(i,m) - cqxx(i,m)
  150           continue
          endif
c
c
c      # modify G fluxes for transverse propagation
c      --------------------------------------------
c
c
c     # split the left-going flux difference into down-going and up-going:
      call rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &          q1d,q1d,aux1,aux2,aux3,
     &          1,amdq,bmasdq,bpasdq)
c
c     # modify flux below and above by B^- A^- Delta q and  B^+ A^- Delta q:
      do 160 m=1,meqn
          do 160 i = 1, mx+1
               gupdate = 0.5d0*dtdx1d(i-1) * bmasdq(i,m)
               gaddm(i-1,m,1) = gaddm(i-1,m,1) - gupdate
               gaddp(i-1,m,1) = gaddp(i-1,m,1) - gupdate
c
               gupdate = 0.5d0*dtdx1d(i-1) * bpasdq(i,m)
               gaddm(i-1,m,2) = gaddm(i-1,m,2) - gupdate
               gaddp(i-1,m,2) = gaddp(i-1,m,2) - gupdate
  160          continue
c
c     # split the right-going flux difference into down-going and up-going:
      call rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &          q1d,q1d,aux1,aux2,aux3,
     &          2,apdq,bmasdq,bpasdq)
c
c     # modify flux below and above by B^- A^+ Delta q and  B^+ A^+ Delta q:
      do 180 m=1,meqn
          do 180 i = 1, mx+1
               gupdate = 0.5d0*dtdx1d(i-1) * bmasdq(i,m)
               gaddm(i,m,1) = gaddm(i,m,1) - gupdate
               gaddp(i,m,1) = gaddp(i,m,1) - gupdate
c
               gupdate = 0.5d0*dtdx1d(i-1) * bpasdq(i,m)
               gaddm(i,m,2) = gaddm(i,m,2) - gupdate
               gaddp(i,m,2) = gaddp(i,m,2) - gupdate
  180          continue
c
  999 continue
      return
      end
