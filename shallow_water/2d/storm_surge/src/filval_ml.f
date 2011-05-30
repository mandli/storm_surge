c
c ------------------------------------------------------------------
c
      subroutine filval(val,mitot,mjtot,hx,hy,lev,time,
     1                  valc,auxc,mic,mjc,
     2                  xleft,xright,ybot,ytop,nvar,
     3                  mptr,ilo,ihi,jlo,jhi,aux,naux,locflip,
     4                  sp_over_h)

      use multilayer_module, only: rho,eta,layers
      use geoclaw_module

      implicit double precision (a-h,o-z)

      include "call.i"

      dimension   val(mitot,mjtot,nvar), valc(mic,mjc,nvar)
      dimension   aux(mitot,mjtot,naux), auxc(mic,mjc,naux)

      double precision coarseval(3,layers)
      logical fineflag(nvar)
      
      double precision :: finemass(2)
      double precision :: slopex(2),slopey(2)

c
c :::::::::::::::::::::::::::::: FILVAL ::::::::::::::::::::::::::
c
c create and fill coarser (lev-1) patch with one extra coarse cell all
c around, plus the ghost cells . will interpolate from this patch to grid mptr
c without needing special boundary code.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     # indext into eta array for surface values:
c       iaddeta(i,j) = loceta + i-1 + mic*(j-1)

      print *,"in filval"

      levc    = lev - 1
      lratiox = intratx(levc)
      lratioy = intraty(levc)
      hxcrse  = hx*lratiox
      hycrse  = hy*lratioy
      xl      = xleft  - hxcrse
      xr      = xright + hxcrse
      yb      = ybot   - hycrse
      yt      = ytop   + hycrse
c
c     set integer indices for coarser patch enlarged by 1 cell
c     (can stick out of domain). proper nesting will insure this one
c     call is sufficient.
      iclo   = ilo/lratiox - 1
      jclo   = jlo/lratioy - 1
      ichi   = (ihi+1)/lratiox - 1 + 1
      jchi   = (jhi+1)/lratioy - 1 + 1
      ng     = 0

c    :::  mcapa  is the capacity function index

      if (naux .eq. 0) then
c     if (mcapa .eq. 0) then
        if (xperdom .or. yperdom .or. spheredom) then
          call preintcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,
     &                    locflip)
        else
          call intcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,1,1)
        endif
      else  ! intersect grids and copy all (soln and aux)
        if (xperdom .or. yperdom .or. spheredom) then
          call preicall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,
     &                  jchi,levc,locflip)
        else
          call icall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi,
     &               levc,1,1)
        endif
      endif
      call bc2amr(valc,auxc,mic,mjc,nvar,naux,
     1            hxcrse,hycrse,levc,time,
     2            xl,xr,yb,yt,
     3            xlower,ylower,xupper,yupper,
     4            xperdom,yperdom,spheredom)

c-----------------------------
c     # for shallow water over topography,
c     # in coarse cells convert from h,
c     # to eta,  before interpolating:
c-----------------------------
      toldry = drytolerance
      
      debug = 42
      open(unit=debug,file='fort.filval', status="unknown")
      
c     #prepare slopes - use min-mod limiters
      do j=2, mjc-1
      do i=2, mic-1
         fineflag(1) = .false.
         
*        !interpolate eta to find depth---------------------------------------
         do ii=-1,1
            if (valc(i+ii,j,4) > toldry) then
                coarseval(2+ii,2) = valc(i+ii,j,4) / rho(2) 
     &               + auxc(i+ii,j,1)
                coarseval(2+ii,1) = valc(i+ii,j,1) / rho(1) 
     &               + coarseval(2+ii,2)
            else
                coarseval(2+ii,2) = eta(2)
                if (valc(i+ii,j,1) > toldry) then
                    coarseval(2+ii,1) = valc(i+ii,j,1) / rho(1) 
     &                   + auxc(i+ii,j,1)
                else
                    coarseval(2+ii,1) = eta(1)
                endif
            endif
         enddo
         
         do m=1,layers
             s1p=coarseval(3,m)-coarseval(2,m)
             s1m=coarseval(2,m)-coarseval(1,m)
             slopex(m)=dmin1(dabs(s1p),dabs(s1m))*dsign(1.d0,
     &                          coarseval(3,m)-coarseval(1,m))
             if (s1m*s1p.le.0.d0) slopex(m) = 0.d0
         enddo

         do jj=-1,1
             if (valc(i,j+jj,4) > toldry) then
                 coarseval(2+jj,2) = valc(i,j+jj,4) / rho(2) 
     &                                  + auxc(i,j+jj,1)
                 coarseval(2+jj,1) = valc(i,j+jj,1) / rho(1) 
     &                                  + coarseval(2+jj,2)
             else
                 coarseval(2+jj,2) = eta(2)
                 if (valc(i,j+jj,1) > toldry) then
                     coarseval(2+jj,1) = valc(i,j+jj,1) 
     &                                      + auxc(i,j+jj,1)
                 else
                     coarseval(2+jj,1) = eta(1)
                 endif
            endif
         enddo

         do m=1,layers
             s1p=coarseval(3,m)-coarseval(2,m)
             s1m=coarseval(2,m)-coarseval(1,m)
             slopey(m)=dmin1(dabs(s1p),dabs(s1m))*dsign(1.d0,
     &                      coarseval(3,m)-coarseval(1,m))
             if (s1m*s1p.le.0.d0) slopey(m)=0.d0
         enddo
         
         
c       !interp. from coarse cells to fine grid to find depth
        finemass = 0.d0
        do ico = 1,lratiox
            do jco = 1,lratioy
                 yoff = (float(jco) - .5)/lratioy - .5
                 xoff = (float(ico) - .5)/lratiox - .5
                 jfine = (j-2)*lratioy + nghost + jco
                 ifine = (i-2)*lratiox + nghost + ico
                 
                 ! Interpolation with eta
                 val(ifine,jfine,4) = coarseval(2,1) + xoff*slopex(2)
     &                                       + yoff*slopey(2)
                 val(ifine,jfine,1) = coarseval(2,1) + xoff*slopex(1) 
     &                                       + yoff*slopey(1)
                 
                 ! Turning val back into h_2
                 val(ifine,jfine,4) = max(0.d0,val(ifine,jfine,4) 
     &                                       - aux(ifine,jfine,1))
                 if (val(ifine,jfine,4) > toldry) then
                     val(ifine,jfine,1) = val(ifine,jfine,1) 
     &                      - val(ifine,jfine,4) - aux(ifine,jfine,1)
                 else
                     fineflag(4) = .true.
                     val(ifine,jfine,5:6) = 0.d0
                     val(ifine,jfine,1) = max(0.d0,val(ifine,jfine,1)
     &                                      -aux(ifine,jfine,1))
                 endif
                 
                 finemass(1) = finemass(1) + val(ifine,jfine,1)
                 finemass(2) = finemass(2) + val(ifine,jfine,4)
                 
                 if (val(ifine,jfine,1) <= toldry) then
                     fineflag(1) = .true.
                     val(ifine,jfine,2) = 0.d0
                     val(ifine,jfine,3) = 0.d0
                 endif
                 
                 ! Reset to actual mass of layer
C                  val(ifine,jfine,1) = rho(1) * val(ifine,jfine,1)
C                  val(ifine,jfine,4) = rho(2) * val(ifine,jfine,4)
            enddo
        enddo
            
* !------determine momentum----------------------------------
*        !finemass is the total mass in all new fine grid cells
*        !all fine mass has been determined for this coarse grid cell
*        !if all fine cells are dry, momentum has already been set
        do m=1,layers
        if (finemass(m) >= toldry) then
            do ivar = 3*(m-1)+2,3*m
               fineflag(ivar)=.false.
               
               s1p=valc(i+1,j,ivar)/rho(m)-valc(i,j,ivar)/rho(m)
               s1m=valc(i,j,ivar)/rho(m)-valc(i-1,j,ivar)/rho(m)
               slopex(m)=dmin1(dabs(s1p),dabs(s1m))*dsign(1.d0,
     &            valc(i+1,j,ivar)/rho(m)-valc(i-1,j,ivar)/rho(m))
               if (s1m*s1p.le.0.d0) slopex(m)=0.d0
               
               s1p=valc(i,j+1,ivar)/rho(m)-valc(i,j,ivar)/rho(m)
               s1m=valc(i,j,ivar)/rho(m)-valc(i,j-1,ivar)/rho(m)
               slopey(m)=dmin1(dabs(s1p),dabs(s1m))*dsign(1.d0,
     &            valc(i,j+1,ivar)/rho(m)-valc(i,j-1,ivar)/rho(m))
               if (s1m*s1p <= 0.d0) slopey(m)=0.d0
               
               if (valc(i,j,3*(m-1)+1)/rho(m) > toldry) then
                  velmax = valc(i,j,ivar)/valc(i,j,3*(m-1)+1)
                  velmin = valc(i,j,ivar)/valc(i,j,3*(m-1)+1)
               else
                  velmax = 0.d0
                  velmin = 0.d0
               endif
               
               do ii = -1,1,2
                  if (valc(i+ii,j,3*(m-1)+1)/rho(m) > toldry) then
                     vel = valc(i+ii,j,ivar)/valc(i+ii,j,3*(m-1)+1)
                     velmax = max(vel,velmax)
                     velmin = min(vel,velmin)
                  endif
                  if (valc(i,j+ii,3*(m-1)+1)/rho(m) > toldry) then
                     vel = valc(i,j,ivar)/valc(i,j+ii,3*(m-1)+1)
                     velmax = max(vel,velmax)
                     velmin = min(vel,velmin)
                  endif
               enddo

*              !try to set momentum
               do ico = 1,lratiox
                  if (fineflag(3*(m-1)+1).or.fineflag(ivar)) exit
                  do jco = 1,lratioy
                     jfine = (j-2)*lratioy + nghost + jco
                     ifine = (i-2)*lratiox + nghost + ico
                     yoff = (float(jco) - .5)/lratioy - .5
                     xoff = (float(ico) - .5)/lratiox - .5
                     hvf = valc(i,j,ivar) / rho(m) + xoff*slopex(m) 
     &                          + yoff * slopey(m)
                     print *,fineflag(3*(m-1)+1),fineflag(ivar)
                     print *,m,3*(m-1)+1,finemass(m)
                     print *,val(ifine,jfine,3*(m-1)+1)
                     print *,hvf
                     vf = hvf/val(ifine,jfine,3*(m-1)+1)
                     if (vf.gt.velmax.or.vf.lt.velmin) then
                        fineflag(ivar)=.true.
                        exit
                     else
                        val(ifine,jfine,ivar) = hvf
                     endif
                     enddo
                  enddo

*              !momentum is set to preserve old momentum or not violate
*              !generating new extrema in velocities
               if (fineflag(1).or.fineflag(ivar)) then !more mass now, conserve momentum
                  area = dble(lratiox*lratioy)
                  dividemass = max(finemass(m),valc(i,j,1))
                  Vnew = area*valc(i,j,ivar)/dividemass

                  do ico = 1,lratiox
                     do jco = 1,lratioy
                        jfine = (j-2)*lratioy + nghost + jco
                        ifine = (i-2)*lratiox + nghost + ico
                        val(ifine,jfine,ivar) = Vnew*val(ifine,jfine,1)
                        enddo
                     enddo
                  endif

               enddo
            endif
            enddo

            ! Multiply by the approrpriate density to get the real
            ! conserved quantities
            do ml=1,layers
                do m=1,3
                    ivar = 3*(ml-1)+m
                    write(debug,*),ivar,val(ifine,jfine,ivar)
                    val(ifine,jfine,ivar) = rho(ml)
     &                       * val(ifine,jfine,ivar)
                enddo
            enddo

         enddo !end of coarse loop
         enddo !end of coarse loop

c
c      if (mcapa .ne. 0) then
c        call fixcapaq(val,aux,mitot,mjtot,valc,auxc,mic,mjc,
c     &                nvar,naux,levc)
c      endif
c
c  overwrite interpolated values with fine grid values, if available.
c
      call intcopy(val,mitot,mjtot,nvar,ilo-nghost,ihi+nghost,
     &             jlo-nghost,jhi+nghost,lev,1,1)

c
c    scan for max wave speed on newly created grid. this will be used to set appropriate
c    time step and appropriate refinement in time. For this app not nec to refine by same
c    amount as refinement in space since refinement at shores where h is shallow has lower
c    speeds.
c
      if (varRefTime) then   ! keep consistent with setgrd_geo and qinit_geo
         sp_over_h = get_max_speed(val,mitot,mjtot,nvar,aux,naux,nghost,
     &                           hx,hy)
      endif

 99   return
      end
