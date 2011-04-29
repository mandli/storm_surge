c -------------------------------------------------------------------
      subroutine flag2refine(mx,my,mbc,meqn,maux,xlower,ylower,dx,dy,
     &                 t,level,tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)
c -------------------------------------------------------------------

c
c ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
c
c User routine to control flagging of points for refinement.
c
c Specific for GeoClaw for tsunami applications and related problems
c 
c
c The logical function allowflag(x,y,t) is called to 
c check whether further refinement at this level is allowed in this cell
c at this time.  
c
c    q   = grid values including ghost cells (bndry vals at specified
c          time have already been set, so can use ghost cell values too)
c
c  aux   = aux array on this grid patch
c
c amrflags  = array to be flagged with either the value
c             DONTFLAG (no refinement needed)  or
c             DOFLAG   (refinement desired)    
c

c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      use hurricane_module
      use topo_module
      use geoclaw_module
      use dtopo_module

      implicit double precision (a-h, o-z)

      dimension   q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      dimension   aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
      dimension   amrflags(1-mbc:mx+mbc,1-mbc:my+mbc)
      logical     allowflag
      external  allowflag
      logical shoreregion,wave,shoreline
      
      double precision :: R_eye(2)

      include 'regions.i'
      include 'qinit.i'

c     # loop over interior points on this grid:
      dmax_grad_wind = 0.d0

      do 200 j = 1,my
        y = ylower +  (j-0.5d0)*dy
        y1 = ylower + (j-1)*dy
        y2 = ylower + j*dy
        do 100 i = 1,mx
          x = xlower +  (i-0.5d0)*dx
          x1 = xlower +  (i-1)*dx
          x2 = xlower +  i*dx
c         # (i,j) grid cell is [x1,x2] x [y1,y2].

c         # default for each point is not to flag unless some condition
c         # below is satisfied:

          amrflags(i,j) = DONTFLAG
          
c         Refine based solely on distance from the eye of the hurricane
          R_eye = t * hurricane_velocity + R_eye_init
          do m=1,max_R_nest
              if((abs(x-R_eye(1)) < R_refine(m)).and.
     &           (abs(y-R_eye(2)) < R_refine(m)).and.
     &           (level <= m)) then
                  amrflags(i,j) = DOFLAG
              endif
          enddo
          
c         Refine based on wind speed
          wind_speed = sqrt(aux(i,j,4)**2+aux(i,j,5)**2)
          do m=1,max_wind_nest
              if ((wind_speed > wind_refine(m)).and.(level <= m)) then
                  amrflags(i,j) = DOFLAG
                  go to 100
              endif
          enddo

          do 30 m=1,mtopofiles
c           # check to see if refinement is forced in any topo file region:
            if (level .lt. minleveltopo(m) .and.
     &          t.ge.tlowtopo(m) .and. t.le.thitopo(m)) then
              xlow = xlowtopo(m)
              xhi = xhitopo(m)
              ylow = ylowtopo(m)
              yhi = yhitopo(m)
                 if (x2.gt.xlow.and.x1.lt.xhi.and.
     &               y2.gt.ylow.and.y1.lt.yhi) then
                    amrflags(i,j) = DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                    endif
              endif
   30       continue

           do 40 m=1,mregions
c           # check to see if refinement is forced in any other region:
            if (level .lt. minlevelregion(m) .and.
     &          t.ge.tlowregion(m) .and. t.le.thiregion(m)) then
              xlow = xlowregion(m)
              xhi = xhiregion(m)
              ylow = ylowregion(m)
              yhi = yhiregion(m)
                 if (x2.gt.xlow.and.x1.lt.xhi.and.
     &               y2.gt.ylow.and.y1.lt.yhi) then
                    amrflags(i,j) = DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                    endif
              endif
   40       continue

         do m = 1,num_dtopo
c           # check if we're in the dtopo region and need to refine:
c           # force refinement to level minleveldtopo
            t0dt = t0dtopo(m)
            tfdt = tfdtopo(m)
            minlevldt = minleveldtopo(m)
            if (level.lt.minleveldtopo(m).and.
     &              t.le.tfdtopo(m).and. !t.ge.t0dtopo(m).and.
     &              x2.gt.xlowdtopo(m).and.x1.lt.xhidtopo(m).and.
     &              y2.gt.ylowdtopo(m).and.y1.lt.yhidtopo(m)) then
                amrflags(i,j)=DOFLAG
                go to 100 !# flagged, so no need to check anything else
                endif
         enddo


         if (iqinit.gt.0 .and. t.eq.0.d0) then
c           # check if we're in the region where initial perturbation is
c           # specified and need to force refinement:
            if (level.lt.minlevelqinit.and.
     &              x2.gt.xlowqinit.and.x1.lt.xhiqinit.and.
     &              y2.gt.ylowqinit.and.y1.lt.yhiqinit) then
                amrflags(i,j)=DOFLAG
                go to 100 !# flagged, so no need to check anything else
                endif
             endif

c        -----------------------------------------------------------------

c        # refinement not forced, so check if it is allowed, and if so,
c        # check if there is a reason to flag this point:

         if (allowflag(x,y,t,level)) then

             surface = q(i,j,1) + aux(i,j,1)
c
c            # RJL: not sure why this next line is done?  
c            # Need to fix for arb.  sealevel?
c            surface = dsign(surface,q(i,j,1))

c            # DLG: it was a way to prevent refining on dry land...
c            # probably should be changed if we allow arbitrary sealevel
c            # by adding sealevel to surface or something.

c            # determine region type and refinement params
            
             shoreregion = dabs(aux(i,j,1)) .lt. depthdeep
             wave = (dabs(surface-sealevel).gt.wavetolerance.and.
     &                q(i,j,1).gt.drytolerance)
c             #DLG: changing following: didn't work so well for non-lat-lon grids
c              shoretol = depthdeep*(dx*dy)
               shoretol = depthdeep
c

             if (wave) then
c               # the surface is not at sea level here

                if (level.lt.maxleveldeep) then
c                   # in deep water we can refine to this level
                    amrflags(i,j)=DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                    endif

c                if (shoreregion) then
c                  shoreline=.false.
c                 # check if any neighboring cell is dry:
c                  do jj=-1,1
c                   do ii=-1,1
c                    shoreline = shoreline.or.q(i+ii,j+jj,1).le.shoretol
c                   enddo
c                  enddo

c                 shoreline=shoreline.and.q(i,j,1).gt.shoretol
                 shoreline = shoreregion

                 if (shoreline.and.q(i,j,1).gt.drytolerance) then
c                    # following comment is regarding commented nested do loop above.
c                    # this cell is wet and a neighbor is dry ==> shoreline
                     amrflags(i,j)=DOFLAG
                     go to 100 !# flagged, so no need to check anything else
                 endif

c               endif
             endif
          endif
          
c        Refine based on momentum or speed of water
c        Had to remove this form the allow flag block as it checks for t > 0
c        and was not allowing refinement before t = 0, we need this as the 
c        storm surge has ramp up time that may need refinement (KTM 2010-8-4)
         speed = sqrt(q(i,j,2)**2 + q(i,j,3)**2)
         if (.not.momentum_refinement) then
             if (q(i,j,1) > drytolerance) then
                 speed = speed / q(i,j,1)
             endif
         endif
         do m=1,max_speed_nest
             if ((speed > speed_refine(m)).and.(level <= m)) then
                 amrflags(i,j) = DOFLAG
                 go to 100
             endif
         enddo

 100     continue  !# end loop on i
 200    continue   !# end loop on j
 
C       print *,dmax_grad_wind

      return
      end
