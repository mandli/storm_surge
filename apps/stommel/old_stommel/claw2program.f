      program driver
c
c  Generic driver routine for claw2
c
c  Author: Randall J. LeVeque
c  Version of March, 1999 --  CLAWPACK Version 4.0
c
c
      implicit double precision (a-h,o-z)

c     # set parameters for maximum array sizes used in declarations
c     # these must be increased for larger problems.
c
c
      parameter (maxmx =   200)
      parameter (maxmy =   200)
      parameter (mwork =  50000)

      parameter (mbc = 2)
      parameter (meqn = 1)
      parameter (mwaves = 1)
      parameter (maux = 3)

      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
      dimension mthlim(mwaves)
      dimension work(mwork)
c
      call claw2ez(maxmx,maxmy,meqn,mwaves,mbc,maux,mwork,mthlim,
     &		 q,work,aux)

      stop 
      end
c
c
c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise
c
       implicit double precision (a-h,o-z)
       dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       common/cdisc/ x0,y0,alf,beta,r0,idisc
       common /stommel/ dlam,bb,dd,capf,rr,ff,gamma,alfs,
     &		       capa,capb,pp,qq,dgbp2

c
       do 20 i=1,mx
	  xi = xlower + (i-0.5d0)*dx
          do 20 j=1,my
	     yj = ylower + (j-0.5d0)*dy
c	     if (xi.lt.0.5d0) then
             xbar = alf*(xi-x0) + beta*(yj-y0)
             ybar = -beta*(xi-x0) + alf*(yj-y0)
	     q(i,j,1) = dcos(xbar*2d-8)*dcos(ybar*2d-8)
	     if (q(i,j,1) .gt. 0.d0) then
 		     q(i,j,1) = 1.d0
 		   else
 		     q(i,j,1) = 0.d0
 		   endif


c	     if (xi.lt.0.5d9 .and. xi.gt.0.1d9 .and. yj.gt.0.1d9 .and.
c    &		 yj.lt.0.3d9) then
c		     q(i,j,1) = 1.d0
c		   else
c		     q(i,j,1) = 0.d0
c		   endif
  20         continue
       return
       end
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann-solver for the advection equation
c     #    q_t  +  u*q_x + v*q_y = 0
c     # where u and v are a given velocity field.
c
c       -----------------------------------------------------------
c     # In advective form, with interface velocities specified by auxl.
c       -----------------------------------------------------------
c
c     # (u,v) may depend of (x,y,t), and the common block comxyt
c     # is used to determine the value of t and either x or y,
c     # depending on what direction the slice is along.
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
      implicit double precision (a-h,o-z)
c
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension apdq(1-mbc:maxm+mbc, meqn)
      dimension amdq(1-mbc:maxm+mbc, meqn)
      dimension auxl(1-mbc:maxm+mbc, *)
      dimension auxr(1-mbc:maxm+mbc, *)
c
c
c     # Set wave, speed, and flux differences:
c     ------------------------------------------
c
      do 30 i = 2-mbc, mx+mbc
         wave(i,1,1) = ql(i,1) - qr(i-1,1)
	 s(i,1) = auxl(i,ixy)
c        # The flux difference df = s*wave  all goes in the downwind direction:
	 amdq(i,1) = dmin1(auxl(i,ixy), 0.d0) * wave(i,1,1)
	 apdq(i,1) = dmax1(auxl(i,ixy), 0.d0) * wave(i,1,1)
   30    continue
c
      return
      end
c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &			imp,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit double precision(a-h,o-z)
c
c     # Riemann solver in the transverse direction for the advection equation.
c
      dimension     ql(1-mbc:maxm+mbc, meqn)
      dimension     qr(1-mbc:maxm+mbc, meqn)
      dimension   asdq(1-mbc:maxm+mbc, meqn)
      dimension bmasdq(1-mbc:maxm+mbc, meqn)
      dimension bpasdq(1-mbc:maxm+mbc, meqn)
      dimension   aux1(1-mbc:maxm+mbc, 2)
      dimension   aux2(1-mbc:maxm+mbc, 2)
      dimension   aux3(1-mbc:maxm+mbc, 2)
c
c
      kv = 3-ixy  !#  = 1 if ixy=2  or  = 2 if ixy=1
      do 10 i=2-mbc,mx+mbc
	 i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
	 bmasdq(i,1) = dmin1(aux2(i1,kv), 0.d0) * asdq(i,1)
	 bpasdq(i,1) = dmax1(aux3(i1,kv), 0.d0) * asdq(i,1)
   10    continue
c
      return
      end
      subroutine setprob
      implicit double precision (a-h,o-z)
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /stommel/ dlam,bb,dd,capf,rr,ff,gamma,alfs,
     &		       capa,capb,pp,qq,dgbp2

c     open(unit=7,file='setprob.data',status='old',form='formatted')
c
c
c     # wall:
      x0 = 0.d0
      y0 = 0.d0
      alf = 1.d0
      beta = 0.d0
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
c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # Set the aux array:
c
c     # aux(i,j,1) = average normal velocity along left edge
c     # aux(i,j,2) = average normal velocity along bottom edge
c     # aux(i,j,3) = capacity of cell: fraction of cell in flow domain
c     #              with some fixes for small triangular cells
c
      implicit real*8(a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 3)
c 
c     # local storage  assumes mbc = 2, maxmx and maxmy < maxm2 - 2
      parameter (maxm2 = 202)
      dimension frleft(-1:maxm2,-1:maxm2), frbot(-1:maxm2,-1:maxm2)
c
c
      do 20 j=1-mbc,my+mbc
         do 20 i=1-mbc,mx+mbc
c
c           # coordinates of lower-left corner and adjacent corners:
            x0 = xlower + (i-1)*dx
            y0 = ylower + (j-1)*dy
            x1 = x0 + dx
            y1 = y0 + dy
c
c           # compute average normal velocity across each edge by
c           # differencing stream function:
c
            aux(i,j,1) = (stream(x0,y0) - stream(x0,y1)) / dy
            aux(i,j,2) = (stream(x1,y0) - stream(x0,y0)) / dx
	    write(17,1701) i,j,aux(i,j,1),aux(i,j,2)
 1701       format(2i4,2d16.6)
c
c           # compute fraction of cell area, left edge, bottom edge,
c           # interior to physical domain
c
            call cartfr(x0,y0,dx,dy,frarea,frleft(i,j),frbot(i,j))
c           if (frarea.lt.0.1d0) frarea = 0.1d0
            aux(i,j,3) = frarea
c
   20       continue

c
c     # triangular cell fix-up:
c
      do 30 j=1-mbc,my+mbc-1
         do 30 i=1-mbc,mx+mbc-1
            hmax = dmax1(frleft(i,j),frbot(i,j),
     &                   frleft(i+1,j),frbot(i,j+1))
            if (hmax .gt. 0.d0) then
                aux(i,j,3) = aux(i,j,3) / hmax
              endif
   30       continue
c
      do  j=1-mbc,my+mbc
         do  i=1-mbc,mx+mbc
	    aux(i,j,3) = dmax1(aux(i,j,3), 1d-3)
	    write(12,*) i,j,aux(i,j,3)
	    enddo
	    enddo
c
      return
      end

c
c
c
c
c     =================================================
      subroutine cartfr(xlow,ylow,dx,dy,wl,hll,hbl)
c     =================================================
      implicit double precision (a-h,o-z)
      external fss
      logical fl(5),alll,allr
      dimension x(10),y(10),xx(5),yy(5)
      common/fsscorn/ xc0,yc0,xc1,yc1
c
c     # compute wl, fraction of cell that lies in left state.
c     # For initial data with two states ql and qr separated by a
c     # discontinuity. The curve along which the discontinuity lies is
c     # specified by the function fdisc, which should return a value that
c     # is negative on the side where ql lies and positive on the qr side.
c
c     # xlow,ylow is the coordinate of the lower left corner of the cell.
c     # dx, dy are grid spacing in x and y.
c
      xx(1) = xlow
      xx(2) = xlow
      xx(3) = xlow+dx
      xx(4) = xlow+dx
      xx(5) = xx(1)
      yy(1) = ylow
      yy(2) = ylow+dy
      yy(3) = ylow+dy
      yy(4) = ylow
      yy(5) = yy(1)
      alll = .true.
      allr = .true.
c
      do 20 i=1,4
         fl(i) = fdisc(xx(i),yy(i)) .lt. 0.d0
         alll = alll .and. fl(i)
         allr = allr .and. (.not. fl(i))
   20    continue
      fl(5) = fl(1)
c
      if (alll) then
         wl = 1.d0
         hll = 1.d0
         hbl = 1.d0
         return
         endif
      if (allr) then
         wl = 0.d0

         hll = 0.d0
         hbl = 0.d0
         return
         endif
c
      iv = 0
      do 40 i=1,4
          if (fl(i)) then
               iv = iv+1
               x(iv) = xx(i)
               y(iv) = yy(i)
               endif
          if (fl(i).neqv.fl(i+1)) then
               iv = iv+1
               xc0 = xx(i)
               yc0 = yy(i)
               xc1 = xx(i+1)
               yc1 = yy(i+1)
               ss = zeroin(0.d0, 1.d0, fss, 1d-8)
c              write(27,*) 'xc,yc,ss:',xc0,yc0,xc1,yc1,ss
               x(iv) = xx(i) + ss*(xx(i+1)-xx(i))
               y(iv) = yy(i) + ss*(yy(i+1)-yy(i))
               endif
c
          if (i.eq.1) then
             if (fl(i) .and. fl(i+1)) then
                 hll = 1.d0
                else if ((.not.fl(i)) .and. (.not.fl(i+1))) then
                 hll = 0.d0
                else if (fl(i)) then
                 hll = ss
                else
                 hll = 1.d0 - ss
                endif
              endif
c
          if (i.eq.4) then
             if (fl(i) .and. fl(i+1)) then
                 hbl = 1.d0
                else if ((.not.fl(i)) .and. (.not.fl(i+1))) then
                 hbl = 0.d0
                else if (fl(i)) then
                 hbl = ss
                else
                 hbl = 1.d0 - ss
                endif
              endif
c
   40     continue
c
c     # compute area:
c
      if (iv.eq.0) then
         wl = 0.d0
         return
         endif
c
      x(iv+1) = x(1)
      y(iv+1) = y(1)
      area = 0.d0
      do 50 i=1,iv
         area = area + .5d0*(y(i)+y(i+1))*(x(i+1)-x(i))
c        write(27,*) '  x,y:',x(i),y(i)
   50    continue
c
      wl = area / (dx*dy)
c     write(27,*) 'area,wl:',area,wl
c
      return
      end
c

c
c
c
c
c     =================================================
      double precision function fss(s)
c     =================================================
      implicit double precision (a-h,o-z)
      common/fsscorn/ xc0,yc0,xc1,yc1
c   
c     # compute fdisc at distance s between corners (xc0,yc0) and (xc1,yc1)
c
      x = xc0 + s*(xc1-xc0)
      y = yc0 + s*(yc1-yc0)
      fss = fdisc(x,y)
      return
      end
c
c
c
c     =================================================
      function zeroin(ax,bx,f,tol)                                         
c     =================================================
      implicit double precision (a-h,o-z)
      double precision f
      external f
c                                                                               
c      a zero of the function  f(x)  is computed in the interval ax,bx .        
c      (Standard routine from netlib)
c                                                                               
c  input..                                                                      
c                                                                               
c  ax     left endpoint of initial interval                                     
c  bx     right endpoint of initial interval                                    
c  f      function subprogram which evaluates f(x) for any x in                 
c         the interval  ax,bx                                                   
c  tol    desired length of the interval of uncertainty of the                  
c         final result ( .ge. 0.0)                                              
c                                                                               
c                                                                               
c  output..                                                                     
c                                                                               
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx             
c                                                                               
c                                                                               
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs        
c  without  a  check.  zeroin  returns a zero  x  in the given interval         
c  ax,bx  to within a tolerance  4*macheps*dabs(x) + tol, where macheps          
c  is the relative machine precision.                                           
c      this function subprogram is a slightly  modified  translation  of        
c  the algol 60 procedure  zero  given in  richard brent, algorithms for        
c  minimization without derivatives, prentice - hall, inc. (1973).              
c                                                                               
c                                                                               
c                                                                               
c  compute eps, the relative machine precision                                  
c                                                                               
      eps = 1.0                                                                 
   10 eps = eps/2.0                                                             
      tol1 = 1.0 + eps                                                          
      if (tol1 .gt. 1.0) go to 10                                               
c                                                                               
c initialization                                                                
c                                                                               
      a = ax                                                                    
      b = bx                                                                    
      fa = f(a)                                                                 
      fb = f(b)                                                                 
c                                                                               
c begin step                                                                    
c                                                                               
   20 c = a                                                                     
      fc = fa                                                                   
      d = b - a                                                                 
      e = d                                                                     
   30 if (dabs(fc) .ge. dabs(fb)) go to 40                                        
      a = b                                                                     
      b = c                                                                     
      c = a                                                                     
      fa = fb                                                                   
      fb = fc                                                                   
      fc = fa                                                                   
c                                                                               
c convergence test                                                              
c                                                                               
   40 tol1 = 2.0*eps*dabs(b) + 0.5*tol                                           
      xm = .5*(c - b)                                                           
      if (dabs(xm) .le. tol1) go to 90                                           
      if (fb .eq. 0.0) go to 90                                                 
c                                                                               
c is bisection necessary                                                        
c                                                                               
      if (dabs(e) .lt. tol1) go to 70                                            
      if (dabs(fa) .le. dabs(fb)) go to 70                                        
c                                                                               
c is quadratic interpolation possible                                           
c                                                                               
      if (a .ne. c) go to 50                                                    
c                                                                               
c linear interpolation                                                          
c                                                                               
      s = fb/fa                                                                 
      p = 2.0*xm*s                                                              
      q = 1.0 - s                                                               
      go to 60                                                                  
c                                                                               
c inverse quadratic interpolation                                               
c                                                                               
   50 q = fa/fc                                                                 
      r = fb/fc                                                                 
      s = fb/fa                                                                 
      p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))                              
      q = (q - 1.0)*(r - 1.0)*(s - 1.0)                                         
c                                                                               
c adjust signs                                                                  
c                                                                               
   60 if (p .gt. 0.0) q = -q                                                    
      p = dabs(p)                                                                
c                                                                               
c is interpolation acceptable                                                   
c                                                                               
      if ((2.0*p) .ge. (3.0*xm*q - dabs(tol1*q))) go to 70                       
      if (p .ge. dabs(0.5*e*q)) go to 70                                         
      e = d                                                                     
      d = p/q                                                                   
      go to 80                                                                  
c                                                                               
c bisection                                                                     
c                                                                               
   70 d = xm                                                                    
      e = d                                                                     
c                                                                               
c complete step                                                                 
c                                                                               
   80 a = b                                                                     
      fa = fb                                                                   
      if (dabs(d) .gt. tol1) b = b + d                                           
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)                              
      fb = f(b)                                                                 
      if ((fb*(fc/dabs(fc))) .gt. 0.0) go to 20                                  
      go to 30                                                                  
c                                                                               
c done                                                                          
c                                                                               
   90 zeroin = b                                                                
      return                                                                    
      end                                                                       
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
c     # into right state
c
      fdisc = (x-x0)*alf + (y-y0)*beta
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

c
c
c
      subroutine claw2ez(maxmx,maxmy,meqn,mwaves,mbc,maux,mwork,mthlim,
     &                   q,work,aux)
c
c     An easy-to-use clawpack driver routine for simple applications
c
c     Author: Randall J. LeVeque
c     Version of August, 1999 --  CLAWPACK Version 4.0
c
      implicit double precision (a-h,o-z)
      external bc2,rpn2,rpt2,src2,b4step2

      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension work(mwork)
      dimension mthlim(mwaves)
c
      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(4)
      dimension tout(100)
      logical outt0
c
      open(55,file='claw2ez.data',status='old',form='formatted')
      open(10,file='fort.info',status='unknown',form='formatted')
      open(11,file='fort.nplot',status='unknown',form='formatted')
c
c
c     # Read the input in standard form from input.claw:
c     # See input.doc for description of each input variable

c     domain variables
      read(55,*) mx
      read(55,*) my

c     i/o variables
      read(55,*) nout
      read(55,*) outstyle
        if (outstyle.eq.1) read(55,*) tfinal
        if (outstyle.eq.2) read(55,*) (tout(i), i=1,nout)
        if (outstyle.eq.3) read(55,*) nsteps


c     timestepping variables
      read(55,*) dtv(1)
      read(55,*) dtv(2)
      read(55,*) cflv(1)
      read(55,*) cflv(2)
      read(55,*) nv(1)
c
	 

c     # input parameters for clawpack routines
      read(55,*) method(1)
      read(55,*) method(2)
      read(55,*) method(3)
      read(55,*) method(4)
      read(55,*) method(5)
      read(55,*) method(6)  
      read(55,*) method(7) 

      read(55,*) meqn1
      read(55,*) mwaves1
      read(55,*) (mthlim(mw), mw=1,mwaves)

      read(55,*) t0
      read(55,*) xlower
      read(55,*) xupper
      read(55,*) ylower
      read(55,*) yupper
c
      read(55,*) mbc1
      read(55,*) mthbc(1)
      read(55,*) mthbc(2)
      read(55,*) mthbc(3)
      read(55,*) mthbc(4)

      if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &    (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
	 write(6,*) '*** ERROR ***  periodic boundary conditions require'
	 write(6,*) '  mthbc(1) and mthbc(2) BOTH be set to 2'
	 stop 
	 endif

      if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or.
     &    (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then
	 write(6,*) '*** ERROR ***  periodic boundary conditions require'
	 write(6,*) '  mthbc(3) and mthbc(4) BOTH be set to 2'
	 stop 
	 endif

c     # These values were passed in, but check for consistency:
c
      if (method(7) .ne. maux) then
         write(6,*) '*** ERROR ***  method(7) should equal maux'
         stop
         endif
      if (meqn1 .ne. meqn) then
	 write(6,*) '*** ERROR ***  meqn set wrong in input or driver'
	 stop
	 endif
      if (mwaves1 .ne. mwaves) then
	 write(6,*) '*** ERROR ***  mwaves set wrong in input or driver'
	 stop
	 endif
      if (mbc1 .ne. mbc) then
	 write(6,*) '*** ERROR ***  mbc set wrong in input or driver'
	 stop
	 endif
c
c     # check that enough storage has been allocated:
c
      if (method(5).lt.2) then
          narray = 1   !# only need one qwork array
        else
          narray = 2   !# need two qwork arrays for Strang splitting
        endif

      maxm = max0(maxmx, maxmy)
      mwork1 = (maxm+2*mbc)*(10*meqn + mwaves + meqn*mwaves 
     &                      + 3*maux + 2) 
     &          + narray * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn   
c
      if (mx.gt.maxmx .or. my.gt.maxmy .or. mwork.lt.mwork1) then
c        # insufficient storage
         maxmx1 = max0(mx,maxmx)
         maxmy1 = max0(my,maxmy)
         maxm1 = max0(maxmx1,maxmy1)

         mwork1 = (maxm1+2*mbc)*(10*meqn + mwaves + meqn*mwaves
     &                      + 3*maux + 2)
     &          + narray * (maxmx1 + 2*mbc) * (maxmy1 + 2*mbc) * meqn

         write(6,*) ' '
         write(6,*) '*** ERROR *** Insufficient storage allocated'
         write(6,*) 'Recompile after increasing values in driver.f:'
         write(6,611) maxmx1
         write(6,612) maxmy1
         write(6,613) mwork1
 611     format(/,'parameter (maxmx = ',i5,')')
 612     format('parameter (maxmy = ',i5,')')
 613     format('parameter (mwork = ',i7,')',/)
         stop
         endif

c
c
      write(6,*) 'running...'
      write(6,*) ' '
c
c     # grid spacing
      dx = (xupper - xlower) / float(mx)
      dy = (yupper - ylower) / float(my)
c


c     # time increments between outputing solution:
      if (outstyle .eq. 1) then
         dtout = (tfinal - t0)/float(nout)
         endif
c
      write(11,1101) nout  
      write(11,1101) 1    
c
 1101 format(i5)
c        
c     # call user's routine setprob to set any specific parameters
c     # or other initialization required.
c
      call setprob
c        
c     # set aux array:
c
      if (maux .gt. 0)  then
         call setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &               maux,aux)
         endif
c
c     # set initial conditions:
c
      call qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &		dx,dy,q,maux,aux)
c
      outt0 = .true.
      if (outt0) then
c        # output initial data
         call out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &          q,t0,0)
         endif

c
c     ----------
c     Main loop:
c     ----------
c
      tend = t0
      do 100 n=1,nout
         tstart = tend
         if (outstyle .eq. 2) then
              tend = tout(n)
            else
              tend = tstart + dtout
            endif
c
         call claw2(maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &           q,aux,xlower,ylower,dx,dy,tstart,tend,dtv,
     &           cflv,nv,method,mthlim,mthbc,
     &           work,mwork,info,bc2,rpn2,rpt2,src2,b4step2)
c
c        # check to see if an error occured:
c        if (info .ne. 0) then
c           write(6,*) '*** ERROR in claw2 ***  info =',info
c           go to 999
c           endif
c
         dtv(1) = dtv(5)  !# use final dt as starting value on next call
c
c        # output solution at this time
c        ------------------------------
c
         call out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &             q,tend,n)
c
c
c        # write out information about this call to claw:
c
         write(6,601) n,tend
  601    format('CLAW2EZ: Frame ',i4,
     &           ' matlab plot files done at time t =',
     &           d12.4,/)
c
         write(10,1010) tend,info,dtv(3),dtv(4),dtv(5),
     &           cflv(3),cflv(4),nv(2)
 1010    format('tend =',d15.4,/,
     &       'info =',i5,/,'smallest dt =',d15.4,/,'largest dt =',
     &       d15.4,/,'last dt =',d15.4,/,'largest cfl =',
     &         d15.4,/,'last cfl =',d15.4,/,'steps taken =',i4,/)
c
  100    continue
c
  999 continue
c
      return 
      end

c
c
c     =====================================================
      subroutine bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &		     dx,dy,q,maux,aux,t,mthbc)
c     =====================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (top), 4 (bottom):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c     ------------------------------------------------
c
c     # Extend the data from the interior cells (1:mx, 1:my)
c     # to the ghost cells outside the region:
c     #   (i, 1-jbc)   for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (i, my+jbc)  for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (1-ibc, j)   for ibc = 1,mbc,  j = 1-mbc, my+mbc
c     #   (mx+ibc, j)  for ibc = 1,mbc,  j = 1-mbc, my+mbc
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
      dimension mthbc(4)

c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc2'
      stop
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do 115 m=1,meqn
         do 115 ibc=1,mbc
            do 115 j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(1,j,m)
  115       continue
      go to 199

  120 continue
c     # periodic:  
      do 125 m=1,meqn
         do 125 ibc=1,mbc
            do 125 j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(mx+1-ibc,j,m)
  125       continue
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 m=1,meqn
         do 135 ibc=1,mbc
            do 135 j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(ibc,j,m)
  135       continue
c     # negate the normal velocity:
      do 136 ibc=1,mbc
         do 136 j = 1-mbc, my+mbc
            q(1-ibc,j,2) = -q(ibc,j,2)
  136    continue
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do 215 m=1,meqn
         do 215 ibc=1,mbc
            do 215 j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(mx,j,m)
  215       continue
      go to 299

  220 continue
c     # periodic:  
      do 225 m=1,meqn
         do 225 ibc=1,mbc
            do 225 j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(ibc,j,m)
  225       continue
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 m=1,meqn
         do 235 ibc=1,mbc
            do 235 j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(mx+1-ibc,j,m)
  235       continue
c     # negate the normal velocity:
      do 236 ibc=1,mbc
         do 236 j = 1-mbc, my+mbc
            q(mx+ibc,j,2) = -q(mx+1-ibc,j,2)
  236    continue
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
      stop
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 m=1,meqn
         do 315 jbc=1,mbc
            do 315 i = 1-mbc, mx+mbc
               q(i,1-jbc,m) = q(i,1,m)
  315       continue
      go to 399

  320 continue
c     # periodic:  
      do 325 m=1,meqn
         do 325 jbc=1,mbc
            do 325 i = 1-mbc, mx+mbc
               q(i,1-jbc,m) = q(i,my+1-jbc,m)
  325       continue
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 m=1,meqn
         do 335 jbc=1,mbc
            do 335 i = 1-mbc, mx+mbc
               q(i,1-jbc,m) = q(i,jbc,m)
  335       continue
c     # negate the normal velocity:
      do 336 jbc=1,mbc
         do 336 i = 1-mbc, mx+mbc
            q(i,1-jbc,3) = -q(i,jbc,3)
  336    continue
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      go to (400,410,420,430) mthbc(4)+1
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
      stop
      go to 499

  410 continue
c     # zero-order extrapolation:
      do 415 m=1,meqn
         do 415 jbc=1,mbc
            do 415 i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,my,m)
  415       continue
      go to 499

  420 continue
c     # periodic:  
      do 425 m=1,meqn
         do 425 jbc=1,mbc
            do 425 i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,jbc,m)
  425       continue
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 m=1,meqn
         do 435 jbc=1,mbc
            do 435 i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,my+1-jbc,m)
  435       continue
c     # negate the normal velocity:
      do 436 jbc=1,mbc
         do 436 i = 1-mbc, mx+mbc
            q(i,my+jbc,3) = -q(i,my+1-jbc,3)
  436    continue
      go to 499

  499 continue

      return
      end
c
c
c =========================================================
      subroutine out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     & 		       dx,dy,q,t,iframe)
c =========================================================
c
c     # Output the results for a general system of conservation laws
c     # in 2 dimensions
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      character*10 fname1, fname2
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw2.m
c     # The same format is used by the amrclaw package.  
c     # Here it's adapted to output just the single grid.
c
c     # first create the file name and open file
c
         fname1 = 'fort.qxxxx'
         fname2 = 'fort.txxxx'
         nstp = iframe
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            fname1(ipos:ipos) = char(ichar('0') + idigit)
            fname2(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue

         open(unit=50,file=fname1,status='unknown',form='formatted')
         open(unit=60,file=fname2,status='unknown',form='formatted')

c
c     # the following parameters are used in amrclaw where there are
c     # multiple grids.  Here they are all set to 1:
      ngrids = 1
      mptr = 1
      level = 1

      write(50,1001) mptr,level,mx,my
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')

      write(50,1002) xlower,ylower,dx,dy
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy',/)
c
      do 20 j=1,my
        do 10 i=1,mx
          do m=1,meqn
c            # exponents with more than 2 digits cause problems reading
c            # into matlab... reset tiny values to zero:
             if (dabs(q(i,j,m)) .lt. 1d-99) q(i,j,m) = 0.d0
             enddo
c
          write(50,1005) (q(i,j,m), m=1,meqn)
 1005     format(4e16.8)
c
 10       continue
        write(50,*) ' '
 20     continue
      write(50,*) ' '

      write(60,1000) t,meqn,ngrids
 1000 format(e18.8,'    time', /, 
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,/)
c

      close(unit=50)
      close(unit=60)

      return
      end
c
c
c
c     ==============================================================
      subroutine claw2(maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &           q,aux,xlower,ylower,dx,dy,tstart,tend,dtv,
     &           cflv,nv,method,mthlim,mthbc,
     &           work,mwork,info,bc2,rpn2,rpt2,src2,b4step2)
c     ==============================================================
c
c
c
c  Solves a hyperbolic system of conservation laws in two space dimensions
c  of the general form
c  
c     capa * q_t + A q_x + B q_y = psi
c
c  The "capacity function" capa(x,y) and source term psi are optional 
c  (see below).
c
c  For a more complete description see the documentation at
c      http://www.amath.washington.edu/~claw
c
c  Sample driver programs and user-supplied subroutines are available.
c  See the the directories claw/clawpack/2d/example* for some examples, and
c  codes in claw/applications for more extensive examples.
c
c  --------------------------------------------------------
c
c  The user must supply the following subroutines:
c
c    bc2, rpn2, rpt2,        subroutines specifying the boundary conditions
c                            and Riemann solvers.
c
c    b4step2            The routine b4step2 is called each time step and
c                       can be supplied by the user in order to perform
c                       other operations that are necessary every time
c                       step.  For example, if the variables stored in
c                       the aux arrays are time-dependent then these
c                       values can be set.   
c
c  In addition, if the equation contains source terms psi, then the user
c  must provide:
c
c    src2               subroutine that solves capa * q_t = psi
c                       over a single time step.
c
c  These routines must be declared EXTERNAL in the main program.
c  For description of the calling sequences, see below.
c
c  Dummy routines b4step1.f and src1.f are available in
c       claw/clawpack/1d/lib
c
c
c
c  Description of parameters...
c  ----------------------------
c
c    maxmx is the maximum number of interior grid points in x, 
c          and is used in declaration of the array q
c
c    maxmy is the maximum number of interior grid points in y, 
c          and is used in declaration of the array q
c
c    meqn is the number of equations in the system of
c         conservation laws.
c
c    mwaves is the number of waves that result from the
c           solution of each Riemann problem.  Often mwaves = meqn but
c           for some problems these may be different, e.g. for the Euler
c           equations meqn = 4 but mwaves = 3 since there are only 3
c           distinct wave speeds.
c
c    mbc is the number of "ghost cells" that must be added on to each
c       side of the domain to handle boundary conditions.  The cells
c       actually in the physical domain are labelled from 1 to mx in x and
c       from 1 to my in y.  The arrays are dimensioned actually indexed
c       from 1-mbc to mx+mbc and from 1-mbc to my+mbc.
c       For the methods currently implemented, mbc = 2 should be used.
c       If the user implements another method that has a larger stencil and
c       hence requires more ghost cells, a larger value of mbc could be used.
c       q is extended from the physical domain to the ghost cells by the
c       user-supplied routine bc2.
c
c    mx is the number of grid cells in the x-direction, in the
c       physical domain.  In addition there are mbc grid cells
c       along each edge of the grid that are used for boundary
c       conditions.
c       Must have mx .le. maxmx
c 
c    my is the number of grid cells in the y-direction, in the
c       physical domain.  In addition there are mbc grid cells
c       along each edge of the grid that are used for boundary
c       conditions.
c       Must have my .le. maxmy
c 
c    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn) 
c        On input:  initial data at time tstart.
c        On output: final solution at time tend.
c        q(i,j,m) = value of mth component in the (i,j) cell.
c        Values within the physical domain are in q(i,j,m) 
c                for i = 1,2,...,mx   and j = 1,2,...,my.
c        mbc extra cells on each end are needed for boundary conditions
c        as specified in the routine bc2.
c
c    aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
c        Array of auxiliary variables that are used in specifying the problem.
c        If method(7) = 0 then there are no auxiliary variables and aux
c                         can be a dummy variable.
c        If method(7) = maux > 0 then there are maux auxiliary variables
c                         and aux must be dimensioned as above.
c
c        Capacity functions are one particular form of auxiliary variable.
c        These arise in some applications, e.g. the
c        determinant of the Jacobian if a mapped grid is used, or a density
c        or porosity function in some advection problems.  
c        See Clawpack Note # 5 for examples.
c
c        If method(6) = 0 then there is no capacity function.
c        If method(6) = mcapa > 0  then there is a capacity function and 
c            capa(i,j), the "capacity" of the (i,j) cell, is assumed to be 
c            stored in aux(i,j,mcapa).
c            In this case we require method(7).ge.mcapa.
c
c    dx = grid spacing in x.  
c         (for a computation in ax <= x <= bx,  set dx = (bx-ax)/mx.)
c
c    dy = grid spacing in y.  
c         (for a computation in ay <= y <= by,  set dy = (by-ay)/my.)
c
c    tstart = initial time.
c
c    tend = Desired final time (on input).
c         = Actual time reached (on output).
c
c    dtv(1:5) = array of values related to the time step:
c               (Note: method(1)=1 indicates variable size time steps)
c         dtv(1) = value of dt to be used in all steps if method(1) = 0
c                = value of dt to use in first step if method(1) = 1
c         dtv(2) = unused if method(1) = 0.
c                = maximum dt allowed if method(1) = 1.
c         dtv(3) = smallest dt used (on output)
c         dtv(4) = largest dt used (on output)
c         dtv(5) = dt used in last step (on output)
c
c    cflv(1:4) = array of values related to Courant number:
c         cflv(1) = maximum Courant number to be allowed.  
c                   With variable time steps the step is retracted and a 
c                   smaller step taken if the Courant
c                   number is larger than this value.  
c                   With fixed time steps the routine aborts.
c                   Usually cflv(1)=1.0 should work 
c                   (or cflv(1)=0.5 if method(3)=0).
c         cflv(2) = unused if method(1) = 0.
c                 = desired Courant number if method(1) = 1.
c                   Should be somewhat less than cflv(1), e.g. 0.9
c         cflv(3) = largest Courant number observed (on output).
c         cflv(4) = Courant number in last step (on output).
c
c    nv(1:2) = array of values related to the number of time steps:
c         nv(1) = unused if method(1) = 0
c               = maximum number of time steps allowed if method(1) = 1
c         nv(2) = number of time steps taken (on output).
c
c    method(1:7) = array of values specifying the numerical method to use
c                  and also indicating whether source terms, capacity
c                  function, auxiliary variables are present in the equation.
c
c         method(1) = 0 if fixed size time steps are to be taken.
c                       In this case, dt = dtv(1) in all steps.
c                   = 1 if variable time steps are to be used.
c                       In this case, dt = dtv(1) in the first step and
c                       thereafter the value cflv(2) is used to choose the
c                       next time step based on the maximum wave speed seen
c                       in the previous step.  Note that since this value
c                       comes from the previous step, the Courant number will
c                       not in general be exactly equal to the desired value
c                       If the actual Courant number in the next step is
c                       greater than cflv(1), then this step is redone with a 
c                       smaller dt.
c
c         method(2) = 1 if only first order increment waves are to be used.
c                   = 2 if second order correction terms are to be added, with
c                       a flux limiter as specified by mthlim.  
c                   = 3 if "third order" correction terms are to be added,
c                       based on my paper "A high-resolution
c                       conservative algorithm for advection in
c                       incompressible flow".
c                       This is still experimental and is currently 
c                       recommended only for problems
c                       with smooth solutions, using no limiter (mthlim = 0)
c
c
c         method(3) = 0 if no transverse propagation is to be applied.
c                       Increment and perhaps correction waves are propagated
c                       normal to the interface.
c                   = 1 if transverse propagation of increment waves 
c                       (but not correction waves, if any) is to be applied.
c                   = 2 if transverse propagation of correction waves is also
c                       to be included.  
c
c                   = -1 if dimensional splitting is to be used instead
c                        of the multi-dimensional wave-propagation.  The
c                        Godunov splitting is used which consists of
c                        sweeping first in x and then in y, with a step of
c                        length dt in each.  The routine bc2 is called
c                        before either sweep to set boundary data, and in
c                        the x-sweep goes over the rows of ghost cells too
c                        so that proper boundary conditions should be set
c                        for the y-sweeps by this process.  Dimensional
c                        splitting is somewhat faster than the unsplit
c                        method and works as well for many (though not all)
c                        problems.
c
c                   = -2 if dimensional splitting is to be used with the
c                        Strang splitting, consisting of 
c                           sweep in x over time dt/2
c                           sweep in y over time dt
c                           sweep in x over time dt/2
c                        This is not recommended because it is slower than
c                        the Godunov splitting and does not appear to be
c                        appreciably better.  Moreover, the boundary
c                        conditions will not be properly set for the final
c                        x-sweep.  (The code could be modified to achieve
c                        this by sweeping over more ghost cells.)
c
c         method(4) = 0 to suppress printing
c                   = 1 to print dt and Courant number every time step
c
c         method(5) = 0 if there is no source term psi.  In this case
c                       the subroutine src2 is never called so a dummy
c                       parameter can be given.
c                   = 1 if there is a source term.  In this case 
c                       the subroutine src2 must be provided and a 
c                       fractional step method is used.
c                       In each time step the following sequence is followed:
c                            call bc to extend data to ghost cells
c                            call step2 to advance hyperbolic eqn by dt
c                            call src2 to advance source terms by dt
c                   = 2 if there is a source term and Strang splitting is to
c                       be used instead of the Godunov splitting above.
c                       In each time step the following sequence is followed:
c                            call bc to extend data to ghost cells
c                            call src2 to advance source terms by dt/2
c                            call step2 to advance hyperbolic equation by dt
c                            call src2 to advance source terms by dt/2
c                       For most problems 1 is recommended rather than 2
c                       since it is less expensive and works essentially as 
c                       well on most problems.  
c                           
c
c         method(6) = 0 if there is no capacity function capa.  
c                   = mcapa > 0 if there is a capacity function.  In this case 
c                       aux(i,j,mcapa) is the capacity of cell (i,j) and you
c                       must also specify method(7) .ge. mcapa and set aux.
c
c         method(7) = 0 if there is no aux array used.
c                   = maux > 0  if there are maux auxiliary variables.
c
c         The recommended choice of methods for most problems is 
c            method(1) = 1,  method(2) = 2,  method(3) = 2.
c
c
c    mthlim(1:mwaves) = array of values specifying the flux limiter to be used
c                     in each wave family mw.  Often the same value will be used
c                     for each value of mw, but in some cases it may be
c                     desirable to use different limiters.  For example,
c                     for the Euler equations the superbee limiter might be
c                     used for the contact discontinuity (mw=2) while another
c                     limiter is used for the nonlinear waves.  Several limiters
c                     are built in and others can be added by modifying the
c                     subroutine philim.
c
c        mthlim(mw) = 0 for no limiter
c                   = 1 for minmod
c                   = 2 for superbee
c                   = 3 for van Leer
c                   = 4 for monotonized centered
c
c    mthbc(1:4) = array of values specifying what boundary conditions should
c                 be used at each edge of the domain, if the standard
c                 bc2.f routine is used.  Passed to bc2.
c
c    work(mwork) = double precision work array of length at least mwork
c
c    mwork = length of work array.  Must be at least
c               N * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn   
c               + (max(mx,my) + 2*mbc) * (10*meqn + mwaves + meqn*mwaves 
c                                          + 3*maux + 2) 
c            where N = 1 if method(5)<2  (no source term or Godunov splitting)
c                  N = 2 if method(5)=2  (source term with Strang splitting)
c            If mwork is too small then the program returns with info = 4
c            and also prints the required value of mwork to unit 6.
c
c            
c    info = output value yielding error information:
c         = 0 if normal return.
c         = 1 if mx.gt.maxmx  or   my.gt.maxmy  or  mbc.lt.2
c         = 2 if method(1)=0 and dt doesn't divide (tend - tstart).
c         = 3 if method(1)=1 and cflv(2) > cflv(1).
c         = 4 if mwork is too small.
c         = 5 if method(6) > method(7)
c         = 11 if the code attempted to take too many time steps, n > nv(1).
c              This could only happen if method(1) = 1 (variable time steps).
c         = 12 if the method(1)=0 and the Courant number is greater than 1
c              in some time step.
c
c           Note: if info.ne.0, then tend is reset to the value of t actually
c           reached and q contains the value of the solution at this time.
c
c    User-supplied subroutines
c    -------------------------
c
c    bc2 = subroutine that specifies the boundary conditions.  
c         This subroutine should extend the values of q from cells
c         (1:mx, 1:my) to the mbc ghost cells along each edge of the domain.
c
c
c    rpn2 = user-supplied subroutine that implements the Riemann solver
c           along a one-dimensional slice of data.
c
c          The form of this subroutine is
c  -------------------------------------------------
c     subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
c    &                  auxl,auxr,wave,s,amdq,apdq)
c
c     implicit double precision (a-h,o-z)
c     dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
c     dimension    s(1-mbc:maxm+mbc, mwaves)
c     dimension   ql(1-mbc:maxm+mbc, meqn)
c     dimension   qr(1-mbc:maxm+mbc, meqn)
c     dimension auxl(1-mbc:maxm+mbc, *)
c     dimension auxr(1-mbc:maxm+mbc, *)
c     dimension amdq(1-mbc:maxm+mbc, meqn)
c     dimension apdq(1-mbc:maxm+mbc, meqn)
c  -------------------------------------------------
c
c         On input, ql contains the state vector at the left edge of each cell
c                   qr contains the state vector at the right edge of each cell
c                 auxl contains auxiliary values at the left edge of each cell
c                 auxr contains auxiliary values at the right edge of each cell
c
c         This data is along a slice in the x-direction if ixy=1
c                                    or the y-direction if ixy=2.
c
c         Note that the i'th Riemann problem has left state qr(i-1,:)
c                                            and right state ql(i,:)
c         In the standard clawpack routines, this Riemann solver is 
c         called with ql=qr=q along this slice.  More flexibility is allowed
c         in case the user wishes to implement another solution method
c         that requires left and rate states at each interface.

c         If method(7)=maux > 0 then the auxiliary variables along this slice
c         are passed in using auxl and auxr.  Again, in the standard routines
c         auxl=auxr is just the values of aux along this slice.

c          On output, 
c             wave(i,m,mw) is the mth component of the jump across
c                              wave number mw in the ith Riemann problem.
c             s(i,mw) is the wave speed of wave number mw in the
c                              ith Riemann problem.
c             amdq(i,m) is the m'th component of the left-going flux difference.
c             apdq(i,m) is the m'th component of the right-going flux difference.
c           It is assumed that each wave consists of a jump discontinuity
c           propagating at a single speed, as results, for example, from a
c           Roe approximate Riemann solver.  An entropy fix can be included
c           into the specification of amdq and apdq.
c
c
c    rpt2 = user-supplied subroutine that implements the splitting of
c           a flux difference asdq into waves in the transverse direction.
c           The form of this subroutine is
c  -------------------------------------------------
c     subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,aux1,aux2,aux3,
c                     imp,asdq,bmasdq,bpasdq)
c
c     implicit double precision (a-h,o-z)
c     dimension     ql(1-mbc:maxm+mbc, meqn)
c     dimension     qr(1-mbc:maxm+mbc, meqn)
c     dimension   aux1(1-mbc:maxm+mbc, maux)
c     dimension   aux2(1-mbc:maxm+mbc, maux)
c     dimension   aux3(1-mbc:maxm+mbc, maux)
c     dimension   asdq(1-mbc:maxm+mbc, meqn)
c     dimension bmasdq(1-mbc:maxm+mbc, meqn)
c     dimension bpasdq(1-mbc:maxm+mbc, meqn)
c  -------------------------------------------------
c          On input, 
c              ql,qr is the data along some one-dimensional slice, as in rpn2
c                   This slice is in the x-direction
c                   if ixy=1, or in the y-direction if ixy=2.  
c              aux2 is the auxiliary array (if method(6)=maux>0) along
c                   this slice, say at j=J if ixy=1.
c              aux1 is the auxiliary array along the adjacent slice J-1
c              aux3 is the auxiliary array along the adjacent slice J+1
c          
c              asdq is an array of flux differences (A^* \Delta q).  
c                   asdq(i,:) is the flux difference propagating away from
c                   the interface between cells i-1 and i.
c              imp = 1 if asdq = A^- \Delta q,  the left-going flux difference
c                    2 if asdq = A^+ \Delta q, the right-going flux difference
c          On output, 
c              bmasdq is the down-going portion of the flux difference
c                   determined by solving a Riemann problem in the transverse
c                   direction using asdq as data.  
c              bpasdq is the up-going portion of the flux difference.
c           For example, for a linear system q_t + Aq_x + Bq_y = 0,  
c                   asdq = A^+ dq  or  A^- dq
c                   and this is then split into
c                       bmasdq = B^- asdq   and   bpasdq = B^+ asdq
c
c
c    src2 = user-supplied subroutine that takes one time step on the 
c           source terms alone, solving
c               capa * q_t = psi
c           over time dt.
c
c           If method(5)=0 then the equation does not contain a source
c           term and this routine is never called.  A dummy argument can
c           be used with many compilers, or provide a dummy subroutine that
c           does nothing (such a subroutine can be found in 
c           clawpack/2d/lib/src2.f)
c
c           The form of this subroutine is
c  -------------------------------------------------
c      subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,q,aux,t,dt)
c      implicit double precision (a-h,o-z)
c      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
c      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
c  -------------------------------------------------
c      If method(7)=0  or the auxiliary variables are not needed in this solver,
c      then the latter dimension statement can be omitted, but aux should
c      still appear in the argument list.
c
c      On input, q(i,j,m) contains the data for solving the 
c                source term equation.
c      On output, q(i,j,m) should have been replaced by the solution to
c                 the source term equation after a step of length dt.
c
c
c
c      b4step2 = subroutine that is called from claw2 before each call to
c                step2.  Use to set time-dependent aux arrays or perform
c                other tasks which must be done every time step.
c
c          The form of this subroutine is
c
c  -------------------------------------------------
c      subroutine b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,
c    &		  xlower,ylower,dx,dy,time,dt,maux,aux)
c      implicit double precision (a-h,o-z)
c      dimension   q(1-mbc:maxmx+mbc, meqn)
c      dimension aux(1-mbc:maxmx+mbc, *)
c  -------------------------------------------------
c
c  
c
c =========================================================================
c
c  Copyright 1994 -- 1999 R. J. LeVeque
c
c  This software is made available for research and instructional use only. 
c  You may copy and use this software without charge for these non-commercial
c  purposes, provided that the copyright notice and associated text is
c  reproduced on all copies.  For all other uses (including distribution of
c  modified versions), please contact the author at the address given below. 
c  
c  *** This software is made available "as is" without any assurance that it
c  *** will work for your purposes.  The software may in fact have defects, so
c  *** use the software at your own risk.
c
c  --------------------------------------
c    CLAWPACK Version 4.0,  August, 1999
c    Webpage: http://www.amath.washington.edu/~claw
c  --------------------------------------
c
c    Author:  Randall J. LeVeque
c             Applied Mathematics
c             Box 352420
c             University of Washington, 
c             Seattle, WA 98195-2420
c             rjl@amath.washington.edu
c =========================================================================
c
c
c
c            
c    ======================================================================
c    Beginning of claw2 code
c    ======================================================================
c 
      implicit double precision (a-h,o-z)
      external bc2,rpn2,rpt2,src2,b4step2
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
      dimension work(mwork)
      dimension mthlim(mwaves),method(7),dtv(5),cflv(4),nv(2),mthbc(4)
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
      maxm = max0(maxmx, maxmy)
      info = 0
      t = tstart
      maxn = nv(1)
      dt = dtv(1)   !# initial dt
      cflmax = 0.d0
      dtmin = dt
      dtmax = dt
      nv(2) = 0
      maux = method(7)
c
c     # check for errors in data:
c     ---------------------------
c
      if (mx.gt.maxmx .or. my.gt.maxmy .or. mbc.lt.2) then
         info = 1
         write(6,*) 'CLAW2 ERROR...  check mx,maxmx,my,maxmy,mbc'
         go to 900
         endif
c
      if (method(1) .eq. 0) then
c        # fixed size time steps.  Compute the number of steps:
         maxn = (tend - tstart + 1d-10) / dt
         if (dabs(maxn*dt - (tend-tstart)) .gt. 1d-8) then
c           # dt doesn't divide time interval integer number of times
            info = 2
            write(6,*) 'CLAW2 ERROR... dt does not divide (tend-tstart)'
            go to 900
            endif
         endif
c
      if (method(1).eq.1 .and. cflv(2).gt.cflv(1)) then
         info = 3
         write(6,*) 'CLAW2 ERROR...  cflv(2) > cflv(1)'
         go to 900
         endif
c
      if (method(6).gt.method(7)) then
         info = 5
         write(6,*) 'CLAW2 ERROR...  method(6) > method(7)'
         go to 900
         endif
c
      if (method(5).lt.2) then
          narray = 1   !# only need one qwork array
        else
          narray = 2   !# need two qwork arrays for Strang splitting
        endif
c
      mwork0 = (maxm+2*mbc)*(10*meqn + mwaves + meqn*mwaves 
     &                      + 3*maux + 2) 
     &          + narray * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn   
c
      if (mwork .lt. mwork0) then
         info = 4
         write(6,*) 'CLAW2 ERROR... mwork should be increased to ',
     &               mwork0
         go to 900
         endif
c
c     # partition work array into pieces needed for local storage in 
c     # step2 routine. Find starting index of each piece:
c
      i0qadd = 1
      i0fadd = i0qadd + (maxm+2*mbc)*meqn
      i0gadd = i0fadd + (maxm+2*mbc)*meqn
      i0q1d = i0gadd + 2*(maxm+2*mbc)*meqn 
      i0dtdx1 = i0q1d + (maxm+2*mbc)*meqn  
      i0dtdy1 = i0dtdx1 + (maxm+2*mbc)
      i0qwrk1 = i0dtdy1 + (maxm+2*mbc)
c
      nqwork = (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn  !# size of q array
      if (method(5).lt.2) then
          i0qwrk2 = i0qwrk1  !# qwrk2 points to same storage as qwrk1
        else
          i0qwrk2 = i0qwrk1 + nqwork  !# second qwork array is needed for
                                      !# Strang spliting
        endif
c
      i0aux1 = i0qwrk2 + nqwork
      i0aux2 = i0aux1 + (maxm+2*mbc)*maux
      i0aux3 = i0aux2 + (maxm+2*mbc)*maux
c
      i0next = i0aux3 + (maxm+2*mbc)*maux  !# next free space
      mused = i0next - 1                  !# space already used
      mwork1 = mwork - mused              !# remaining space (passed to step2)
c
c
c
c     -----------
c     # main loop
c     -----------
c
      if (maxn.eq.0 .or. tend.le.tstart) go to 900
      do 100 n=1,maxn
         told = t   !# time at beginning of time step.
         if (told+dt .gt. tend) dt = tend - told
c
   40    continue
c
c        # store dt and t in the common block comxyt in case they are needed
c        # in the Riemann solvers (for variable coefficients)
         tcom = told
         dtcom = dt
         dxcom = dx
         dycom = dy
c
c
c
c        ================================================================
c
c        -------------------------
c        # main steps in algorithm
c        -------------------------
c
c        # extend data from grid to bordering boundary cells:
         call bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &		     dx,dy,q,maux,aux,told,mthbc)
c
c
c
c        # call user-supplied routine which might set aux arrays
c        # for this time step, for example.

         call b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,
     &                xlower,ylower,dx,dy,told,dt,maux,aux)
c
c
c
         if (method(5).eq.2) then
c            # with Strang splitting for source term:
c            # First need to store solution before taking
c            # step on source terms in case we need to redo everything with
c            # a smaller time step if the Courant number is too large in 
c            # subroutine step2.
             call copyq2(maxmx,maxmy,meqn,mbc,mx,my,q,work(i0qwrk2))
c
c            # source terms over a half time step:
             dt2 = dt / 2.d0
             call src2(maxmx,maxmy,meqn,mbc,mx,my,q,aux,told,dt2)
             endif
c
c        # copy q into qwork1.  q is updated in step2 and qwork1 is
c        # preserved to provide data for Riemann problems.
c        # qwork1 can also be used to restart if the Courant number is 
c        # too large, unless Strang splitting is used in which case we
c        # must used the values already stored above before
c        # taking the source term step.
         call copyq2(maxmx,maxmy,meqn,mbc,mx,my,q,work(i0qwrk1))
c
c        # take one step on the conservation law:
c
         if( method(3) .ge. 0 )then
c            # unsplit version
c             
             call step2(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                  work(i0qwrk1),q,aux,
     &                  dx,dy,dt,method,mthlim,cfl,
     &                  work(i0qadd),work(i0fadd),work(i0gadd),
     &                  work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &                  work(i0aux1),work(i0aux2),work(i0aux3),
     &                  work(i0next),mwork1,rpn2,rpt2)
c
         else
c           # dimensional splitting (fractional steps)
c
            call dimsp2(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                  work(i0qwrk1),q,aux,
     &                  dx,dy,dt,method,mthlim,cfl,cflv,
     &                  work(i0qadd),work(i0fadd),work(i0gadd),
     &                  work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &                  work(i0aux1),work(i0aux2),work(i0aux3),
     &                  work(i0next),mwork1,rpn2,rpt2)
c
         endif
c
         t = told + dt
c
         if (method(4).eq.1) then
c            # verbose mode
             write(6,601) n,cfl,dt,t
  601        format('CLAW2... Step',i6,
     &                   '   Courant number =',f6.3,'  dt =',d12.4,
     &                   '  t =',d12.4)
             endif
c
c
c        # check to see if the Courant number was too large:
         if (cfl .le. cflv(1)) then
c             # accept this step
              cflmax = dmax1(cfl,cflmax)
            else
c             # Reject this step.  Reset q to qwork from previous time:
c             # Note that i0qwrk2 points to work space where previous
c             # solution is stored in all cases method(5) = 0,1, or 2. 
              t = told
              call copyq2(maxmx,maxmy,meqn,mbc,mx,my,work(i0qwrk2),q)
c
              if (method(4).eq.1) then
c                # verbose mode
                 write(6,602)
                 endif
  602         format('CLAW2 rejecting step... Courant number too large')
c
              if (method(1).eq.1) then
c                 # if variable dt, go back and take a smaller step.
                  dt = dmin1(dtv(2), dt * cflv(2)/cfl)
                  go to 40
                else
c                 # if fixed dt, give up and return
                  cflmax = dmax1(cfl,cflmax)
                  go to 900
                endif
            endif
c
c
c        # claw2 step is accepted
c        # now apply source terms:
c
         if (method(5).eq.2) then
c            # source terms over a second half time step for Strang splitting:
c            # Note it is not so clear what time t should be used here if
c            # the source terms are time-dependent!
             call src2(maxmx,maxmy,meqn,mbc,mx,my,q,aux,t,dt2)
             endif
c
         if (method(5).eq.1) then
c            # source terms over a full time step:
             call src2(maxmx,maxmy,meqn,mbc,mx,my,q,aux,t,dt)
             endif
c
c        ================================================================
c
c
c
         if (method(1) .eq. 1) then
c           # choose new time step if variable time step
            if (cfl.eq.0.d0) then 
                dt = dtv(2)
              else
                dt = dmin1(dtv(2), dt * cflv(2)/cfl)
              endif
            dtmin = dmin1(dt,dtmin)
            dtmax = dmax1(dt,dtmax)
            endif
c
c
c        # see if we are done:
c
         nv(2) = nv(2) + 1
         if (t .ge. tend) go to 900
  100    continue
c
  900  continue
c 
c      # return information
c
       if (method(1).eq.1 .and. t.lt.tend .and. nv(2) .eq. maxn) then
c         # too many timesteps
          write(6,*) 'CLAW2 ERROR...  too many timesteps'
          info = 11
          endif
       if (method(1).eq.0 .and. cflmax .gt. cflv(1)) then
c         # Courant number too large with fixed dt
          write(6,*) 'CLAW2 ERROR...  Courant number too large'
          info = 12
          endif
       tend = t
       cflv(3) = cflmax
       cflv(4) = cfl
       dtv(3) = dtmin
       dtv(4) = dtmax
       dtv(5) = dt
c
       return 
       end
c
c
c
c
c     ==========================================================
      subroutine step2(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &               qold,qnew,aux,dx,dy,dt,method,mthlim,cfl,
     &               qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2)
c     ==========================================================
c
c     # Take one time step, updating q.
c     # On entry, qold and qnew should be identical and give the
c     #    initial data for this step
c     # On exit, qnew returns values at the end of the time step.
c     #    qold is unchanged.
c    
c     # qadd is used to return increments to q from flux2
c     # fadd and gadd are used to return flux increments from flux2.
c     # See the flux2 documentation for more information.
c
c
      implicit double precision (a-h,o-z)
      dimension qold(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension qnew(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  q1d(1-mbc:maxm+mbc, meqn)
      dimension qadd(1-mbc:maxm+mbc, meqn)
      dimension fadd(1-mbc:maxm+mbc, meqn)
      dimension gadd(1-mbc:maxm+mbc, meqn, 2)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
      dimension aux1(1-mbc:maxm+mbc, *)
      dimension aux2(1-mbc:maxm+mbc, *)
      dimension aux3(1-mbc:maxm+mbc, *)

      dimension dtdx1d(1-mbc:maxmx+mbc)
      dimension dtdy1d(1-mbc:maxmx+mbc)
      dimension method(7),mthlim(mwaves)
      dimension work(mwork)

      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
c
c
c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c
      i0wave = 1
      i0s = i0wave + (maxm+2*mbc)*meqn*mwaves
      i0amdq = i0s + (maxm+2*mbc)*mwaves
      i0apdq = i0amdq + (maxm+2*mbc)*meqn
      i0cqxx = i0apdq + (maxm+2*mbc)*meqn
      i0bmadq = i0cqxx + (maxm+2*mbc)*meqn
      i0bpadq = i0bmadq + (maxm+2*mbc)*meqn
      iused = i0bpadq + (maxm+2*mbc)*meqn - 1
c
      if (iused.gt.mwork) then
c        # This shouldn't happen due to checks in claw2
         write(6,*) '*** not enough work space in step2'
         write(6,*) '*** iused = ', iused, '   mwork =',mwork
         stop 
      endif
c
c
      mcapa = method(6)
      maux = method(7)
      cfl = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
c
      if (mcapa.eq.0) then
c        # no capa array:
         do 5 i=1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
    5    continue
      endif
c
c
c     # perform x-sweeps
c     ==================
c
      do 50 j = 0,my+1
c
c        # copy data along a slice into 1d arrays:
         do 21 m=1,meqn
            do 20 i = 1-mbc, mx+mbc
               q1d(i,m) = qold(i,j,m)
   20       continue
   21    continue
c
         if (mcapa.gt.0)  then
            do 22 i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(i,j,mcapa)
   22       continue
         endif
c
         if (maux .gt. 0)  then
            do 25 ma=1,maux
               do 24 i = 1-mbc, mx+mbc
                  aux1(i,ma) = aux(i,j-1,ma)
                  aux2(i,ma) = aux(i,j  ,ma)
                  aux3(i,ma) = aux(i,j+1,ma)
   24          continue
   25       continue
         endif
c
c     # Store the value of j along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         jcom = j  
c           
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(1,maxm,meqn,mwaves,mbc,mx,
     &            q1d,dtdx1d,aux1,aux2,aux3,method,mthlim,
     &            qadd,fadd,gadd,cfl1d,
     &              work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &              work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2)
         cfl = dmax1(cfl,cfl1d)
c
c        # update qnew by flux differencing.
c        # (rather than maintaining arrays f and g for the total fluxes,
c        # the modifications are used immediately to update qnew
c        # in order to save storage.)
c
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:
            do 31 m=1,meqn
               do 30 i=1,mx
                  qnew(i,j,m) = qnew(i,j,m) + qadd(i,m)
     &                 - dtdx * (fadd(i+1,m) - fadd(i,m))
     &                       - dtdy * (gadd(i,m,2) - gadd(i,m,1))
                  qnew(i,j-1,m) = qnew(i,j-1,m) - dtdy * gadd(i,m,1)
                  qnew(i,j+1,m) = qnew(i,j+1,m) + dtdy * gadd(i,m,2)
   30          continue
   31       continue
c
         else
c
c            # with capa array.  
            do 41 m=1,meqn
               do 40 i=1,mx
                  qnew(i,j,m) = qnew(i,j,m) + qadd(i,m)
     &                 - (dtdx * (fadd(i+1,m) - fadd(i,m))
     &                         +  dtdy * (gadd(i,m,2) - gadd(i,m,1)))
     &                       / aux(i,j,mcapa)
                 qnew(i,j-1,m) = qnew(i,j-1,m) - dtdy * gadd(i,m,1)
     &                       / aux(i,j-1,mcapa)
                 qnew(i,j+1,m) = qnew(i,j+1,m) + dtdy * gadd(i,m,2)
     &                       / aux(i,j+1,mcapa)
   40          continue
   41       continue
         endif
   50 continue
c
c
c
c     # perform y sweeps
c     ==================
c
c
      do 100 i = 0, mx+1
c
c        # copy data along a slice into 1d arrays:
         do 71 m=1,meqn
            do 70 j = 1-mbc, my+mbc
               q1d(j,m) = qold(i,j,m)
   70       continue
   71    continue
c
         if (mcapa.gt.0)  then
            do 72 j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(i,j,mcapa)
   72       continue
         endif
c
         if (maux .gt. 0)  then
            do 75 ma=1,maux
               do 74 j = 1-mbc, my+mbc
                  aux1(j,ma) = aux(i-1,j,ma)
                  aux2(j,ma) = aux(i,  j,ma)
                  aux3(j,ma) = aux(i+1,j,ma)
   74          continue
   75       continue
         endif
c
c
c     # Store the value of i along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         icom = i  
c           
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(2,maxm,meqn,mwaves,mbc,my,
     &            q1d,dtdy1d,aux1,aux2,aux3,method,mthlim,
     &            qadd,fadd,gadd,cfl1d,
     &              work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &              work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2)
c
         cfl = dmax1(cfl,cfl1d)
c
c        # update qnew by flux differencing.
c        # Note that the roles of fadd and gadd are reversed for
c        # the y-sweeps -- fadd is the modification to g-fluxes and
c        # gadd is the modification to f-fluxes to the left and right.
c
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:
            do 81 m=1,meqn
               do 80 j=1,my
                  qnew(i,j,m) = qnew(i,j,m) + (qadd(j,m)
     &                  - dtdy * (fadd(j+1,m) - fadd(j,m))
     &                  - dtdx * (gadd(j,m,2) - gadd(j,m,1)))
                  qnew(i-1,j,m) = qnew(i-1,j,m) - dtdx * gadd(j,m,1)
                  qnew(i+1,j,m) = qnew(i+1,j,m) + dtdx * gadd(j,m,2)
   80          continue
   81       continue
c
         else
c
c            # with capa array.  
            do 91 m=1,meqn
               do 90 j=1,my
                  qnew(i,j,m) = qnew(i,j,m) + qadd(j,m)
     &                  - (dtdy * (fadd(j+1,m) - fadd(j,m))
     &                    + dtdx * (gadd(j,m,2) - gadd(j,m,1)))
     &                       / aux(i,j,mcapa)
                  qnew(i-1,j,m) = qnew(i-1,j,m) - dtdx * gadd(j,m,1)
     &                       / aux(i-1,j,mcapa)
                  qnew(i+1,j,m) = qnew(i+1,j,m) + dtdx * gadd(j,m,2)
     &                       / aux(i+1,j,mcapa)
   90          continue
   91       continue
         endif
  100 continue
c
      return
      end
c     ============================================
      subroutine b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,
     &		  xlower,ylower,dx,dy,time,dt,maux,aux)
c     ============================================
c
c     # called from claw2 before each call to step2.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.

c
c     # dummy routine 
c
c     
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, *)
c
      return
      end
c
c
c
c
c     ==========================================================
      subroutine step2ds(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &               qold,qnew,aux,dx,dy,dt,method,mthlim,cfl,
     &               qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2,ids)
c     ==========================================================
c
c     # Take one time step, updating q.
c     # On entry, qold and qnew should be identical and give the
c     #    initial data for this step
c     # On exit, qnew returns values at the end of the time step.
c     #    qold is unchanged.
c    
c     # qadd is used to return increments to q from flux2
c     # fadd and gadd are used to return flux increments from flux2.
c     # See the flux2 documentation for more information.
c
c
      implicit double precision (a-h,o-z)
      external rpn2,rpt2
      dimension qold(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension qnew(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  q1d(1-mbc:maxm+mbc, meqn)
      dimension qadd(1-mbc:maxm+mbc, meqn)
      dimension fadd(1-mbc:maxm+mbc, meqn)
      dimension gadd(1-mbc:maxm+mbc, meqn, 2)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
      dimension aux1(1-mbc:maxm+mbc, *)
      dimension aux2(1-mbc:maxm+mbc, *)
      dimension aux3(1-mbc:maxm+mbc, *)

      dimension dtdx1d(1-mbc:maxmx+mbc)
      dimension dtdy1d(1-mbc:maxmx+mbc)
      dimension method(7),mthlim(mwaves)
      dimension work(mwork)

      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
c
c
c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c
      i0wave = 1
      i0s = i0wave + (maxm+2*mbc)*meqn*mwaves
      i0amdq = i0s + (maxm+2*mbc)*mwaves
      i0apdq = i0amdq + (maxm+2*mbc)*meqn
      i0cqxx = i0apdq + (maxm+2*mbc)*meqn
      i0bmadq = i0cqxx + (maxm+2*mbc)*meqn
      i0bpadq = i0bmadq + (maxm+2*mbc)*meqn
      iused = i0bpadq + (maxm+2*mbc)*meqn - 1
c
      if (iused.gt.mwork) then
c        # This shouldn't happen due to checks in claw2
         write(6,*) '*** not enough work space in step2'
         write(6,*) '*** iused = ', iused, '   mwork =',mwork
         stop 
      endif
c
c
      mcapa = method(6)
      maux = method(7)
      cfl = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
c
      if (mcapa.eq.0) then
c        # no capa array:
         do 5 i=1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
    5    continue
      endif
c
      if( ids.eq.1 )then
c
c     # perform x-sweeps
c     ==================
c
c     # note that for dimensional splitting we sweep over the rows of
c     # ghosts cells as well as the interior.  This updates the ghost
c     # cell values to the intermediate state as needed in the following 
c     # sweep in the y-direction.
c
      do 50 j = 1-mbc,my+mbc
c
c        # copy data along a slice into 1d arrays:
         do 21 m=1,meqn
            do 20 i = 1-mbc, mx+mbc
               q1d(i,m) = qold(i,j,m)
   20       continue
   21    continue
c
         if (mcapa.gt.0)  then
            do 22 i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(i,j,mcapa)
   22       continue
         endif
c
         if (maux .gt. 0)  then
             do 23 ma=1,maux
               do 23 i = 1-mbc, mx+mbc
                 aux2(i,ma) = aux(i,j  ,ma)
   23          continue
c
             if(j .ne. 1-mbc)then
                do 24 ma=1,maux
                   do 24 i = 1-mbc, mx+mbc
                      aux1(i,ma) = aux(i,j-1,ma)
   24              continue
                endif
c
             if(j .ne. my+mbc)then
                do 25 ma=1,maux
                   do 25 i = 1-mbc, mx+mbc
                      aux3(i,ma) = aux(i,j+1,ma)
   25              continue
                endif
c
             endif
c
c        # Store the value of j along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         jcom = j  
c           
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(1,maxm,meqn,mwaves,mbc,mx,
     &            q1d,dtdx1d,aux1,aux2,aux3,method,mthlim,
     &            qadd,fadd,gadd,cfl1d,
     &              work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &              work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2)
         cfl = dmax1(cfl,cfl1d)
c
c        # update qnew by flux differencing.
c        # (rather than maintaining arrays f and g for the total fluxes,
c        # the modifications are used immediately to update qnew
c        # in order to save storage.)
c
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:
            do 31 m=1,meqn
               do 30 i=1,mx
                  qnew(i,j,m) = qnew(i,j,m) + qadd(i,m)
     &                 - dtdx * (fadd(i+1,m) - fadd(i,m))
   30          continue
   31       continue
c
         else
c
c            # with capa array.  
            do 41 m=1,meqn
               do 40 i=1,mx
                  qnew(i,j,m) = qnew(i,j,m) + qadd(i,m)
     &                        - dtdx * (fadd(i+1,m) - fadd(i,m))
     &                        / aux(i,j,mcapa)
   40          continue
   41       continue
         endif
   50 continue
c
      endif
c
      if( ids.eq.2 )then
c
c     # perform y sweeps
c     ==================
c
c
      do 100 i = 1-mbc, mx+mbc
c
c        # copy data along a slice into 1d arrays:
         do 71 m=1,meqn
            do 70 j = 1-mbc, my+mbc
               q1d(j,m) = qold(i,j,m)
   70       continue
   71    continue
c
         if (mcapa.gt.0)  then
            do 72 j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(i,j,mcapa)
   72       continue
         endif
c
         if (maux .gt. 0)  then
c
             do 73 ma=1,maux
               do 73 j = 1-mbc, my+mbc
                 aux2(j,ma) = aux(i,j,ma)
   73          continue
c
             if(i .ne. 1-mbc)then
                do 74 ma=1,maux
                   do 74 j = 1-mbc, my+mbc
                      aux1(j,ma) = aux(i-1,j,ma)
   74              continue
                endif
c
             if(i .ne. mx+mbc)then
                do 75 ma=1,maux
                   do 75 j = 1-mbc, my+mbc
                      aux3(j,ma) = aux(i+1,j,ma)
   75              continue
                endif
c
             endif            
c
c     # Store the value of i along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         icom = i  
c           
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(2,maxm,meqn,mwaves,mbc,my,
     &            q1d,dtdy1d,aux1,aux2,aux3,method,mthlim,
     &            qadd,fadd,gadd,cfl1d,
     &              work(i0wave),work(i0s),work(i0amdq),work(i0apdq),
     &              work(i0cqxx),work(i0bmadq),work(i0bpadq),rpn2,rpt2)
c
         cfl = dmax1(cfl,cfl1d)
c
c        # update qnew by flux differencing.
c        # Note that the roles of fadd and gadd are reversed for
c        # the y-sweeps -- fadd is the modification to g-fluxes and
c        # gadd is the modification to f-fluxes to the left and right.
c
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:
            do 81 m=1,meqn
               do 80 j=1,my
                  qnew(i,j,m) = qnew(i,j,m) + qadd(j,m)
     &                  - dtdy * (fadd(j+1,m) - fadd(j,m))
   80          continue
   81       continue
c
         else
c
c            # with capa array.  
            do 91 m=1,meqn
               do 90 j=1,my
                  qnew(i,j,m) = qnew(i,j,m) + qadd(j,m)
     &                  - dtdy * (fadd(j+1,m) - fadd(j,m))
     &                        / aux(i,j,mcapa)
   90          continue
   91       continue
         endif
  100 continue
c
      endif
c
c
      return
      end
c
c
c     ==========================================================
      subroutine dimsp2(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                  qold,qnew,aux,dx,dy,dt,method,mthlim,cfl,
     &                  cflv,qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                  aux1,aux2,aux3,work,mwork,rpn2,rpt2)
c     ==========================================================
c
c     # Take one time step, updating q, using dimensional
c     # splitting. Two choices are available:
c     #
c     # method(3) = -1   gives Godunov splitting:
c     #    time step dt in x-direction
c     #    time step dt in y-direction
c
c     # method(3) = -2   gives Strang splitting
c     #    time step dt/2 in x-direction
c     #    time step dt   in y-direction
c     #    time step dt/2 in x-direction
c
c     # Godunov splitting is recommended over Strang splitting normally
c     # since it typically works as well, is faster, and boundary
c     # conditions are handled properly.
c
      implicit double precision (a-h,o-z)
      dimension qold(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension qnew(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  q1d(1-mbc:maxm+mbc, meqn)
      dimension cflv(4)
      dimension qadd(1-mbc:maxm+mbc, meqn)
      dimension fadd(1-mbc:maxm+mbc, meqn)
      dimension gadd(1-mbc:maxm+mbc, meqn, 2)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
      dimension aux1(1-mbc:maxm+mbc, *)
      dimension aux2(1-mbc:maxm+mbc, *)
      dimension aux3(1-mbc:maxm+mbc, *)

      dimension dtdx1d(1-mbc:maxmx+mbc)
      dimension dtdy1d(1-mbc:maxmx+mbc)
      dimension method(7),mthlim(mwaves)
      dimension work(mwork)

c
c     # If method(3) = -1, take a full time step in x.
c     # If method(3) = -2, take a half time step in x.
c
      dt2 = dt/2.d0
c
      if( method(3) .eq. -2 )then
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                 qold,qnew,aux,dx,dy,dt2,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
      else
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                 qold,qnew,aux,dx,dy,dt,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
      endif
c
      if (cflx .gt. cflv(1)) then
c        # Abort if the Courant number was too large in x-sweep
         cfl = cflx
	 return
	 endif
c
c     # Take full step in y-direction
c
      call step2ds(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &             qnew,qnew,aux,dx,dy,dt,method,mthlim,cfly,
     &             qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &             aux1,aux2,aux3,work,mwork,rpn2,rpt2,2)
c
      cfl = dmax1(cflx,cfly)
c
c     # Finally, take a half time step in the x-direction
c     # if Strang splitting is used.  NOTE: boundary conditions may
c     # not be set properly for this sweep.
c
      if( method(3) .eq. -2 )then
         if (cfly .gt. cflv(1)) then
c           # Abort if the Courant number was too large in y-sweep
            cfl = cfly
	    return
	    endif
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                 qnew,qnew,aux,dx,dy,dt2,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
          cfl = dmax1(cfl,cflx)
      endif
c
      return
      end



c
c
c
c
c     =====================================================
      subroutine flux2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                 q1d,dtdx1d,aux1,aux2,aux3,method,mthlim,
     &               qadd,fadd,gadd,cfl1d,wave,s,
     &                 amdq,apdq,cqxx,bmasdq,bpasdq,rpn2,rpt2)
c     =====================================================
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
c     Note that if method(6)=1 then the capa array comes into the second 
c     order correction terms, and is already included in dtdx1d:
c     If ixy = 1 then
c        dtdx1d(i) = dt/dx                 if method(6) = 0
c                  = dt/(dx*capa(i,jcom))  if method(6) = 1
c     If ixy = 2 then
c        dtdx1d(j) = dt/dy                 if method(6) = 0
c                  = dt/(dy*capa(icom,j))  if method(6) = 1
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
c
      implicit double precision (a-h,o-z)
      external rpn2,rpt2
      dimension    q1d(1-mbc:maxm+mbc, meqn)
      dimension   amdq(1-mbc:maxm+mbc, meqn)
      dimension   apdq(1-mbc:maxm+mbc, meqn)
      dimension bmasdq(1-mbc:maxm+mbc, meqn)
      dimension bpasdq(1-mbc:maxm+mbc, meqn)
      dimension   cqxx(1-mbc:maxm+mbc, meqn)
      dimension   qadd(1-mbc:maxm+mbc, meqn)
      dimension   fadd(1-mbc:maxm+mbc, meqn)
      dimension   gadd(1-mbc:maxm+mbc, meqn, 2)
c
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension aux1(1-mbc:maxm+mbc, *)
      dimension aux2(1-mbc:maxm+mbc, *)
      dimension aux3(1-mbc:maxm+mbc, *)
c
      dimension     s(1-mbc:maxm+mbc, mwaves)
      dimension  wave(1-mbc:maxm+mbc, meqn, mwaves)
c
      dimension method(7),mthlim(mwaves)
      logical limit
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
      limit = .false.
      do 5 mw=1,mwaves
         if (mthlim(mw) .gt. 0) limit = .true.
   5  continue
c
c     # initialize flux increments:
c     -----------------------------
c
         do 20 m=1,meqn
            do 10 i = 1-mbc, mx+mbc
               qadd(i,m) = 0.d0
               fadd(i,m) = 0.d0
               gadd(i,m,1) = 0.d0
               gadd(i,m,2) = 0.d0
   10       continue
   20    continue
c
c
c     # solve Riemann problem at each interface and compute Godunov updates
c     ---------------------------------------------------------------------
c
      call rpn2(ixy,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux2,aux2,
     &        wave,s,amdq,apdq)
c
c     # Set qadd for the donor-cell upwind method (Godunov)
      do 41 m=1,meqn
         do 40 i=1,mx+1
            qadd(i,m) = qadd(i,m) - dtdx1d(i)*apdq(i,m)
            qadd(i-1,m) = qadd(i-1,m) - dtdx1d(i-1)*amdq(i,m)
   40    continue
   41 continue
c
c     # compute maximum wave speed for checking Courant number:
      cfl1d = 0.d0
      do 51 mw=1,mwaves
         do 50 i=1,mx+1
c          # if s>0 use dtdx1d(i) to compute CFL,
c          # if s<0 use dtdx1d(i-1) to compute CFL:
            cfl1d = dmax1(cfl1d, dtdx1d(i)*s(i,mw), 
     &				-dtdx1d(i-1)*s(i,mw))
c          # initialize cqxx (used for second order correction terms)
	   cqxx(i,m) = 0.d0
   50    continue
   51 continue
c
      if (method(2).eq.1) go to 130
c
c     # modify F fluxes for second order q_{xx} correction terms:
c     -----------------------------------------------------------
c
c     # apply limiter to waves:
      if (limit) call limiter(maxm,meqn,mwaves,mbc,mx,wave,s,mthlim)
c
      do 121 m=1,meqn
         do 120 i = 1, mx+1
c
c        # For correction terms below, need average of dtdx in cell
c        # i-1 and i.  Compute these and overwrite dtdx1d:
c
            dtdx1d(i-1) = 0.5d0 * (dtdx1d(i-1) + dtdx1d(i))
c
            cqxx(i,m) = 0.d0
            do 119 mw=1,mwaves
c
c              # second order corrections:
               cqxx(i,m) = cqxx(i,m) + dabs(s(i,mw))
     &            * (1.d0 - dabs(s(i,mw))*dtdx1d(i-1)) * wave(i,m,mw)
c
c             ---------------------------------------------------------
c              # third order corrections:
c              # still experimental... works well for smooth solutions
c              # with no limiters but not well with limiters so far.
c              # dtdx not handled properly for method(6)=1
c
               if (method(2).lt.3) go to 119
               if (s(i,mw) .gt. 0.d0) then
                  dq2 = wave(i,m,mw) - wave(i-1,m,mw)
               else
                  dq2 = wave(i+1,m,mw) - wave(i,m,mw)
               endif
               cqxx(i,m) = cqxx(i,m) - s(i,mw)/3.d0 *
     &                     (1.d0 - (s(i,mw)*dtdx1d(i-1))**2) * dq2
c             ---------------------------------------------------------
c
  119       continue
            fadd(i,m) = fadd(i,m) + 0.5d0 * cqxx(i,m)
  120    continue
  121 continue
c
c
  130  continue
c
      if (method(3).le.0) go to 999   !# no transverse propagation
c
      if (method(3).eq.2) then
c         # incorporate cqxx into amdq and apdq so that it is split also.
         do 151 m=1,meqn
            do 150 i = 1, mx+1
               amdq(i,m) = amdq(i,m) + cqxx(i,m)
               apdq(i,m) = apdq(i,m) - cqxx(i,m)
  150       continue
  151    continue
      endif
c
c
c      # modify G fluxes for transverse propagation
c      --------------------------------------------
c
c
c     # split the left-going flux difference into down-going and up-going:
      call rpt2(ixy,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,aux3,
     &        1,amdq,bmasdq,bpasdq)
c
c     # modify flux below and above by B^- A^- Delta q and  B^+ A^- Delta q:
      do 161 m=1,meqn
         do 160 i = 1, mx+1
            gadd(i-1,m,1) = gadd(i-1,m,1) - 
     &                 0.5d0*dtdx1d(i-1) * bmasdq(i,m)
            gadd(i-1,m,2) = gadd(i-1,m,2) -
     &                 0.5d0*dtdx1d(i-1) * bpasdq(i,m)
  160    continue
  161 continue
c
c     # split the right-going flux difference into down-going and up-going:
      call rpt2(ixy,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,aux3,
     &        2,apdq,bmasdq,bpasdq)
c
c     # modify flux below and above by B^- A^+ Delta q and  B^+ A^+ Delta q:
      do 181 m=1,meqn
         do 180 i = 1, mx+1
            gadd(i,m,1) = gadd(i,m,1) - 
     &                0.5d0*dtdx1d(i-1) * bmasdq(i,m)
            gadd(i,m,2) = gadd(i,m,2) - 
     &                0.5d0*dtdx1d(i-1) * bpasdq(i,m)
  180    continue
  181 continue
c
  999 continue
      return
      end
c
c
c =========================================================
      subroutine copyq2(maxmx,maxmy,meqn,mbc,mx,my,q1,q2)
c =========================================================
c
c     # copy the contents of q1 into q2
c
      implicit double precision (a-h,o-z)
      dimension q1(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension q2(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
c
c
      do 10 j = 1-mbc, my+mbc
        do 10 i = 1-mbc, mx+mbc
          do 10 m=1,meqn
             q2(i,j,m) = q1(i,j,m)
   10        continue
      return
      end
c
c
c     =====================================================
      subroutine limiter(maxm,meqn,mwaves,mbc,mx,wave,s,mthlim)
c     =====================================================
c
c     # Apply a limiter to the waves.
c     # The limiter is computed by comparing the 2-norm of each wave with
c     # the projection of the wave from the interface to the left or
c     # right onto the current wave.  For a linear system this would
c     # correspond to comparing the norms of the two waves.  For a 
c     # nonlinear problem the eigenvectors are not colinear and so the 
c     # projection is needed to provide more limiting in the case where the
c     # neighboring wave has large norm but points in a different direction
c     # in phase space.
c
c     # The specific limiter used in each family is determined by the
c     # value of the corresponding element of the array mthlim, as used in
c     # the function philim.
c     # Note that a different limiter may be used in each wave family.
c
c     # dotl and dotr denote the inner product of wave with the wave to
c     # the left or right.  The norm of the projections onto the wave are then
c     # given by dotl/wnorm2 and dotr/wnorm2, where wnorm2 is the 2-norm
c     # of wave.
c
      implicit real*8(a-h,o-z)
      dimension mthlim(mwaves)
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
c
c
      do 50 mw=1,mwaves
         if (mthlim(mw) .eq. 0) go to 50
         dotr = 0.d0
         do 40 i = 0, mx+1
            wnorm2 = 0.d0
            dotl = dotr
            dotr = 0.d0
            do 20 m=1,meqn
               wnorm2 = wnorm2 + wave(i,m,mw)**2
               dotr = dotr + wave(i,m,mw)*wave(i+1,m,mw)
   20          continue
            if (i.eq.0) go to 40
            if (wnorm2.eq.0.d0) go to 40
c
            if (s(i,mw) .gt. 0.d0) then
                wlimitr = philim(wnorm2, dotl, mthlim(mw))
              else
                wlimitr = philim(wnorm2, dotr, mthlim(mw))
              endif
c
            do 30 m=1,meqn
               wave(i,m,mw) = wlimitr * wave(i,m,mw)
   30          continue
   40       continue
   50    continue
c
      return
      end
c
c
c     =====================================================
      double precision function philim(a,b,meth)
c     =====================================================
      implicit real*8(a-h,o-z)
c
c     # Compute a limiter based on wave strengths a and b.
c     # meth determines what limiter is used.
c     # a is assumed to be nonzero.
c
      r = b/a
      go to (10,20,30,40,50) meth

c
   10 continue
c     --------
c     # minmod
c     --------
      philim = dmax1(0.d0, dmin1(1.d0, r))
      return
c
   20 continue
c     ----------
c     # superbee
c     ----------
      philim = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
      return
c
   30 continue
c     ----------
c     # van Leer
c     ----------
      philim = (r + dabs(r)) / (1.d0 + dabs(r))
      return
c
   40 continue
c     ------------------------------
c     # monotinized centered 
c     ------------------------------
      c = (1.d0 + r)/2.d0
      philim = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
      return
c
   50 continue
c     ------------------------------
c     # Beam-Warming
c     ------------------------------
      philim = r

      return
      end
c
c      =======================================================
       subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &		       dx,dy,q,maux,aux,t,dt)
c      =======================================================
c
       implicit double precision (a-h,o-z)
       dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
c
c      # dummy subroutine for use when equation has no source term.
c      # If method(5)=0 then this routine is never called, but its
c      # existence may be required by some compilers.
c
       return
       end
