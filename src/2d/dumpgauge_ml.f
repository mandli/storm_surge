      subroutine dumpgauge(q,aux,xlow,ylow,nvar,mitot,mjtot,mptr)

      use multilayer_module, only: layers, rho
      use geoclaw_module

      implicit double precision (a-h,o-z)

      include "gauges.i"
      include "call.i"

      integer bsearch
      dimension q(mitot,mjtot,nvar), var(maxvar)
      dimension aux(mitot,mjtot,1)
      dimension eta(layers)
      dimension h(layers,4)

c     See if this grid contains any gauges so data can be output
c     may turn out this should be sorted, but for now do linear search
c     
c     Array is sorted according to indices in mbestorder array
c     so do binary search to find start. Could have many same source grids

      if (mgauges == 0) return

      istart = bsearch(mptr) 
      if (istart < 1) return      ! This grid not used

c     This stuff the same for all gauges on this grid
      tgrid = rnode(timemult,mptr)
      level = node(nestlevel,mptr)
      hx    =  hxposs(level)
      hy    =  hyposs(level)

      do ii = istart, mgauges
          i = mbestorder(ii)   ! gauge number
          if (mptr /= mbestsrc(i)) return  ! all done
          
          if (tgrid < t1gauge(i) .or. tgrid > t2gauge(i)) then
c            Don't output at this time for gauge i
             cycle
          endif
          
c     prepare to do linear interp at gauge location to get vars
c     should fancier (limited) interp be done?
          iindex =  int(.5 + (xgauge(i)-xlow)/hx)
          jindex =  int(.5 + (ygauge(i)-ylow)/hy)
          if ((iindex < nghost .or. iindex > mitot-nghost) .or.
     &       (jindex < nghost .or. jindex > mjtot-nghost)) then
             print *,"ERROR in output of Gauge Data ",i,iindex,jindex
          endif
          xcent = xlow + (iindex-.5)*hx
          ycent = ylow + (jindex-.5)*hy
          xoff  = abs((xgauge(i)-xcent)/hx)
          yoff  = abs((ygauge(i)-ycent)/hy)
	      if (xoff < 0 .or. xoff > 1 .or. yoff < 0. .or. yoff > 1) then
	          print *," BIG PROBLEM in DUMPGAUGE", i,xoff,yoff  
    	  endif

c      Modified by RJL 12/31/09 to interpolate only where all four cells are
c      wet, otherwise just take this cell value:

c      Check for dry cells by comparing h to drytol2, which should be smaller
c      than drytolerance to avoid oscillations since when h < drytolerance the
c      velocities are zeroed out which can then lead to increase in h again.

          drytol2 = 0.1d0 * drytolerance

          do m=1,layers
              layer_index = 3*(m-1)
              h(m,1) = q(iindex,jindex,layer_index+1) / rho(m)
              h(m,2) = q(iindex+1,jindex,layer_index+1) / rho(m)
              h(m,3) = q(iindex,jindex+1,layer_index+1) / rho(m)
              h(m,4) = q(iindex+1,jindex+1,layer_index+1) / rho(m)
              
              if ((h(m,1) < drytol2) .or.
     &            (h(m,2) < drytol2) .or.
     &            (h(m,3) < drytol2) .or.
     &            (h(m,4) < drytol2)) then
                  ! One of the cells is dry, so just use value from grid cell
                  ! that contains gauge rather than interpolating
                  
                  icell = int(1.d0 + (xgauge(i) - xlow) / hx)
                  jcell = int(1.d0 + (ygauge(i) - ylow) / hy)
                  do ivar=1,3
                      var(ivar + layer_index) = 
     &                       q(icell,jcell,ivar + layer_index) / rho(m)
                  enddo
                  if (m == layers) then
                      ! This is the bottom layer and we should figure out the
                      ! topography
                      topo = aux(icell,jcell,1)
                  endif
              else
                  ! Linear interpolation between four cells
                  do ivar=1,3
                      var(layer_index + ivar) = (1.d0 - xoff) * 
     &                   (1.d0 - yoff)
     &                 * q(iindex,jindex,layer_index + ivar) / rho(m)
     &                 + xoff*(1.d0 - yoff) 
     &                 * q(iindex+1,jindex,layer_index + ivar) / rho(m)
     &                 + (1.d0 - xoff) * yoff 
     &                 * q(iindex,jindex+1,layer_index + ivar) / rho(m)
     &                 + xoff * yoff 
     &                 * q(iindex+1,jindex+1,layer_index+ivar) / rho(m)
                  enddo
                  if (m == layers) then
                      topo = (1.d0 - xoff) * (1.d0 - yoff) 
     &                        * aux(iindex,jindex,1) 
     &                      + xoff * (1.d0 - yoff) 
     &                        * aux(iindex+1,jindex,1) 
     &                      + (1.d0 - xoff) * yoff 
     &                        * aux(iindex,jindex+1,1) 
     &                      + xoff * yoff 
     &                        * aux(iindex+1,jindex+1,1)
                  endif
              endif
          enddo
              
          ! Calculate eta, this assumes two layers
          eta(2) = h(2,1) + topo
          eta(1) = h(1,1) + eta(2)
              
          write(OUTGAUGEUNIT,100) igauge(i),level,tgrid, 
     &                    (var(j),j=1,3*layers),(eta(j),j=1,layers)
      enddo
      
 100  format(2i5,15e15.7)
 
      end subroutine dumpgauge
c
c --------------------------------------------------------------------
c
      subroutine setbestsrc()
c
c  ## called every time grids change, to set the best source grid
c  ## to find gauge data
c
c  ## lbase is grid level that didn't change but since fine
c  ## grid may have disappeared, still have to look starting
c  ## at coarsest level 1.
c
      implicit double precision (a-h,o-z)

      include "gauges.i"
      include "call.i"
c
c ##  set source grid for each loc from coarsest level to finest.
c ##  that way finest src grid left and old ones overwritten
c ##  this code uses fact that grids do not overlap

c # for debugging, initialize sources to 0 then check that all set
      do i = 1, mgauges
         mbestsrc(i) = 0
      end do

 
      do 20 lev = 1, lfine  
          mptr = lstart(lev)
 5        do 10 i = 1, mgauges
            if ((xgauge(i) .ge. rnode(cornxlo,mptr)) .and.    
     .          (xgauge(i) .le. rnode(cornxhi,mptr)) .and.    
     .          (ygauge(i) .ge. rnode(cornylo,mptr)) .and.  
     .          (ygauge(i) .le. rnode(cornyhi,mptr)) )
     .      mbestsrc(i) = mptr
 10       continue

          mptr = node(levelptr, mptr)
          if (mptr .ne. 0) go to 5
 20   continue


      do i = 1, mgauges
        if (mbestsrc(i) .eq. 0) 
     .      write(6,*)"ERROR in setting grid src for gauge data",i
      end do

c
c     sort the source arrays for easy testing during integration
      call qsorti(mbestorder,mgauges,mbestsrc)

      return
      end
c
c ------------------------------------------------------------------------
c
      integer function bsearch(mptr)

      implicit double precision (a-h,o-z)
      include "gauges.i"

      bsearch = -1           ! signal if not found

      indexlo = 1
      indexhi = mgauges

 5    if (indexhi .lt. indexlo) go to 99
      mid = (indexlo + indexhi)/2
      mbomid = mbestorder(mid)

      if (mptr .gt. mbestsrc(mbestorder(mid))) then
	   indexlo = mid+1
	   go to 5
      else if (mptr .lt. mbestsrc(mbestorder(mid))) then
	   indexhi = mid-1
	   go to 5
      else    ! found the grid. find its first use in the array
	 istart = mid


 10      if (istart .gt. 1) then
            if (mbestsrc(mbestorder(istart-1)) .ne. mptr) go to 90
            istart = istart - 1
            go to 10
          endif

      endif

 90   bsearch = istart

 99   return
      end
	 
