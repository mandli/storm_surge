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
      use multilayer_module
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
      logical shoreregion,shoreline
 
      logical wave
      dimension wave(2)     
      double precision speed
      dimension speed(2)
      integer layer
      
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
          
c ----------------------------------------------------------------------------
c  Hurricane Refinement criteria
c
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
c ----------------------------------------------------------------------------
c  Refinement based on regions
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
C          if (allowflag(x,y,t,level)) then
         if (.true.) then
             
             wave = .false.
c            Check to see if bottom layer is dry
             if (q(i,j,4) / rho(2) > drytolerance) then
                 eta2 = q(i,j,4) / rho(2) + aux(i,j,1)
                 wave(2) = (abs(eta2-eta(2)) > wave_tol(2))
             else
                 eta2 = aux(i,j,1)
                 wave(2) = .false.
             endif
             eta1 = q(i,j,1) / rho(1) + eta2
             
             shoreregion = abs(aux(i,j,1)) < depthdeep
             wave(1) = (abs(eta1-eta(1)) > wave_tol(1)).and.
     &               (q(i,j,1) / rho(1) > drytolerance)
             
C              wave = (((dabs(eta1-sealevel) > wavetolerance).and.
C      &                (q(i,j,4) / rho(2) > drytolerance)).or.
C      &               ((dabs(eta2-internallevel) > wavetolerance).and.
C      &                (q(i,j,1) / rho(1) > drytolerance)))
     
            
C              shoreregion = dabs(aux(i,j,1)) .lt. depthdeep
C              wave = (dabs(surface-sealevel).gt.wavetolerance.and.
C      &                q(i,j,1).gt.drytolerance)
c             #DLG: changing following: didn't work so well for non-lat-lon grids
c              shoretol = depthdeep*(dx*dy)
               shoretol = depthdeep
c

             if (wave(1).or.wave(2)) then
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

c                
             endif
C            Gradient based refinement
             dq = 0.d0
             do ml=1,layers
                 do m=1,3 
                 index = 3*(m-1)
                 dqi = abs(q(i+1,j,index+m) / rho(ml) 
     &                       - q(i-1,j,index+m) / rho(ml))
                 dqj = abs(q(i,j+1,index+m) / rho(ml) 
     &                       - q(i,j-1,index+m) / rho(ml))
                 dq = max(dq,dqi,dqj)
             enddo
             if (dq > tolsp) amrflags(i,j) = DOFLAG
             
          endif
          
c        Refine based on momentum or speed of water
c        Had to remove this form the allow flag block as it checks for t > 0
c        and was not allowing refinement before t = 0, we need this as the 
c        storm surge has ramp up time that may need refinement (KTM 2010-8-4)
         do m=1,layers
             index = 3*(m-1)
             if (q(i,j,index+1) > drytolerance) then
                 speed(m) = sqrt((q(i,j,index+2) / q(i,j,index+1))**2 
     &                         + (q(i,j,index+3) / q(i,j,index+1))**2)
             else
                 speed(m) = 0.d0
             endif
         enddo

         do m=1,max_speed_nest
             do index=1,layers
                 if ((speed(index) > speed_refine(m)).and.
     &               (level <= m)) then
                     amrflags(i,j) = DOFLAG
                     goto 100
                 endif
             enddo
         enddo

 100    continue  !# end loop on i
 200    continue   !# end loop on j
 
C       print *,dmax_grad_wind

      return
      end
