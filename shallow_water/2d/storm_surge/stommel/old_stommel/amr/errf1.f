c
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,rctcrse,mptr,mi2tot,mj2tot,
     2                 mitot,mjtot,rctflg,sperr)
      implicit double precision (a-h,o-z)

      include  "call.i"
 
      dimension  rctfine(mitot,mjtot,nvar)
      dimension  rctcrse(mi2tot,mj2tot,nvar)
      dimension  rctflg(mitot,mjtot,nvar)
      dimension  sperr(mitot,mjtot)
      logical    allowed
       common/cdisc/ x0,y0,alf,beta,r0,idisc
c
c
c ::::::::::::::::::::::::::::: ERRF1 ::::::::::::::::::::::::::::::::
c
c  compare error estimates in rctfine, rctcrse. If exceed tol, flag.
c  We put in a hook which allows us to ignore the error estimates
c  and not refine if we are outside some prescribed region.
c
c  sperr is the spatial component of the error estimate only
c  provides simple way for user to specify extra flagging in errsp
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     # to allow refinement anywhere:
c     allowed(x,y,level) = .true.
c
c     # This function can be changed to allow refinement only on some portion
c     # of the grid.
c     # This is useful if you wish to zoom in on some structure in a 
c     # known location but don't want the same level of refinement elsewhere.  
c     # Points are flagged only if one of the errors is greater than the 
c     # corresponding tolerance.
c
c
      allowed(x,y,level) = 
     &        (level.le.1 .or. (alf*(x-x0) + beta*(y-y0)).le.0.1d8) 
c    &        (level.le.1 .or. x.le.0.4d8) 
c
 
      time  = rnode(timemult, mptr)
      xleft = rnode(cornxlo,mptr)
      levm  = node(nestlevel, mptr)
      hx    = hxposs(levm)
      ybot  = rnode(cornylo,mptr)
      hy    = hyposs(levm)
      dt    = possk(levm)
      numsp = 0
 
      errmax = 0.0d0
      err2   = 0.0d0
c     order  = dt*dfloat(2**(iorder+1) - 2)
      order  = dfloat(2**(iorder+1) - 2)
c
      if (.not. (edebug)) go to 20
         write(outunit,107) mptr
 107     format(//,' coarsened grid values for grid ',i4)
         do 10 jj = nghost+1, mj2tot-nghost
            j = mj2tot + 1 - jj
            write(outunit,101) (rctcrse(i,j,1),
     .				i = nghost+1, mi2tot-nghost)
10       continue
         write(outunit,108) mptr
 108     format(//, ' fine grid values for grid ',i4)
         do 15 jj = nghost+1, mjtot-nghost
            j = mjtot + 1 - jj
            write(outunit,101) (rctfine(i,j,1),i=nghost+1,mitot-nghost)
15       continue
101      format(' ',13f6.3)
c
c zero out the exterior locations so they don't affect err.est.
c
 20   continue
      jfine = nghost+1
      do 35  j = nghost+1, mj2tot-nghost
      yofj  = ybot + (dfloat(jfine) - .5d0)*hy
      ifine = nghost+1
c
      do 30  i  = nghost+1, mi2tot-nghost
          rflag = goodpt
          xofi  = xleft + (dfloat(ifine) - .5d0)*hx
          term1 = rctfine(ifine,jfine,1)
          term2 = rctfine(ifine+1,jfine,1)
          term3 = rctfine(ifine+1,jfine+1,1)
          term4 = rctfine(ifine,jfine+1,1)
c         # divide by (aval*order) for relative error
          aval  = (term1+term2+term3+term4)/4.d0
          est   =  dabs((aval-rctcrse(i,j,1))/ order)
          if (est .gt. errmax) errmax = est
	  err2 = err2 + est*est
c         write(outunit,102) i,j,est
 102      format(' i,j,est ',2i5,e12.5)
c         rctcrse(i,j,2) = est
c
          if (est .ge. tol .and. allowed(xofi,yofj,levm)) then
             rflag  = badpt
          endif 
      rctcrse(i,j,1) = rflag
      ifine = ifine + 2
 30   continue
      jfine = jfine + 2
 35   continue
c
c  transfer flagged points on cell centered coarse grid
c  to cell centered fine grid. count flagged points.
c
c  initialize rctflg to 0.0 (no flags)  before flagging
c
      do 40 j = 1, mjtot
      do 40 i = 1, mitot
 40      rctflg(i,j,1) = goodpt
c
c  print out intermediate flagged rctcrse (for debugging)
c
      if (eprint) then
	 err2 = dsqrt(err2/dfloat((mi2tot-2*nghost)*(mj2tot-2*nghost)))
         write(outunit,103) mptr, levm, errmax, err2
 103     format(' grid ',i4,' level ',i4,
     .          ' max. error = ',e15.7,' err2 = ',e15.7)
	 if (edebug) then
	   write(outunit,*) ' flagged points on coarsened grid ',
     .             	    'for grid ',mptr
	   do 45 jj = nghost+1, mj2tot-nghost
	      j = mj2tot + 1 - jj
	      write(outunit,106) (nint(rctcrse(i,j,1)),
     .                            i=nghost+1,mi2tot-nghost)
106           format(1h ,80i1)
45         continue
         endif
      endif
c
      jfine   = nghost+1
      do 70 j = nghost+1, mj2tot-nghost
      ifine   = nghost+1
      do 60 i = nghost+1, mi2tot-nghost
	 if (rctcrse(i,j,1) .eq. goodpt) go to 55
	    rctflg(ifine,jfine,1)    = badpt
	    rctflg(ifine+1,jfine,1)  = badpt
	    rctflg(ifine,jfine+1,1)  = badpt
	    rctflg(ifine+1,jfine+1,1)= badpt
 55       ifine   = ifine + 2
 60     continue
        jfine   = jfine + 2
 70   continue
c
      if (edebug) then
	 write(outunit,*)" spatial error for grid ",mptr
         do 75 jjfine = nghost+1, mjtot-nghost
	    jfine = mjtot + 1 - jjfine
	    write(outunit,101)(sperr(ifine,jfine),
     .			       ifine=nghost+1,mitot-nghost)
 75      continue
      endif

      do 80 jfine = nghost+1, mjtot-nghost
      yofj  = ybot + (dfloat(jfine) - nghost - .5d0)*hy
      do 80 ifine = nghost+1, mitot-nghost
        xofi  = xleft + (dfloat(ifine) - nghost - .5d0)*hx
        if (sperr(ifine,jfine) .gt. tolsp .and. allowed(xofi,yofj,levm))
     &     then
	       rflag = rctflg(ifine,jfine,1)
               if (rflag .ne. badpt) then
	         rctflg(ifine,jfine,1) = badpt
      	         numsp = numsp + 1
               endif
	  endif
 80   continue

      if (eprint) then
	 write(outunit,118) numsp,mptr
 118     format( i5,' more pts. flagged for spatial error on grid',i4,/)
	if (edebug) then
	  do 56 jj = nghost+1, mjtot-nghost
	   j = mjtot + 1 - jj
	   write(outunit,106)(nint(rctflg(i,j,1)),i=nghost+1,mitot-nghost)
 56       continue
	endif
      endif

      return
      end
