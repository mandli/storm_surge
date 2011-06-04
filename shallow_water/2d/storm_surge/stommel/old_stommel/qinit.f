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

       ic = 2
       go to (10,20) ic
c
   10  continue
c      # checkerboard:

       do 11 i=1,mx
	  xi = xlower + (i-0.5d0)*dx
          do 11 j=1,my
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
  11         continue
      go to 999

  20  continue
c
c      # blob:

       do 21 i=1,mx
	  xi = xlower + (i-0.5d0)*dx
          do 21 j=1,my
	     yj = ylower + (j-0.5d0)*dy
c	     if (xi.lt.0.5d0) then
             xbar = alf*(xi-x0) + beta*(yj-y0)
             ybar = -beta*(xi-x0) + alf*(yj-y0)
 	     if ((xbar-3.d8)**2 + (ybar-2.d8)**2 .lt. 1.d16) then
 		     q(i,j,1) = 1.d0
 		   else
 		     q(i,j,1) = 0.d0
 		   endif
  21         continue
 999   continue
       return
       end
