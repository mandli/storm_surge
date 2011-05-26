c     ============================================
      subroutine b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,
     &            xlower,ylower,dx,dy,t,dt,maux,aux)
c     ============================================
c
c     # called before each call to step
c     # use to set time-dependent aux arrays or perform other tasks.
c
c     This particular routine sets negative values of q(i,j,1) to zero,
c     as well as the corresponding q(i,j,m) for m=1,meqn.
c     This is for problems where q(i,j,1) is a depth.
c     This should occur only because of rounding error.

c     Also calls movetopo if topography might be moving.

      use multilayer_module
      use hurricane_module
      use geoclaw_module
      use topo_module
      use dtopo_module

      implicit double precision (a-h,o-z)

      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)
      dimension vel(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, 2)
      
      integer :: layer,layer_index
      double precision :: h(2),u(2),v(2),g
      logical :: dry_state(2)
      
      g = grav
      
c=====================Parameters===========================================


c     # check for NANs in solution:
      call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,1)

c     # check for h < 0 and reset to zero
c     # check for h < drytolerance
c     # set hu = hv = 0 in all these cells

      do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
            do m=1,layers
              index = 3*(m-1)
              if (q(i,j,index+1) / rho(m) < drytolerance) then
                 q(i,j,index+1) = max(q(i,j,index+1),0.d0)
                 q(i,j,index+2) = 0.d0
                 q(i,j,index+3) = 0.d0
              endif
            enddo
        enddo
      enddo

      write(26,*) 'B4STEP2: t, num_dtopo: ', t,num_dtopo
      do i=1,num_dtopo
          call movetopo(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,t,dt,maux,aux,
     &      dtopowork(i0dtopo(i):i0dtopo(i)+mdtopo(i)-1),
     &      xlowdtopo(i),ylowdtopo(i),xhidtopo(i),yhidtopo(i),
     &      t0dtopo(i),tfdtopo(i),dxdtopo(i),dydtopo(i),dtdtopo(i),
     &      mxdtopo(i),mydtopo(i),mtdtopo(i),mdtopo(i),
     &      minleveldtopo(i),maxleveldtopo(i),topoaltered(i))
      enddo
    
      ! Set wind and pressure aux variables for this grid
      write(26,*) "B4STEP2:  Setting aux array for wind and pressure"
      call hurricane_wind(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,t,
     & aux(:,:,4:5))
      call hurricane_pressure(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,
     & dy,t,aux(:,:,6))
     
      ! Check Richardson number
      if (layers > 1) then
      do i=1,mx
          do j=1,my
              dry_state = .false.
              do layer=1,2
                  m = 3*(layer-1)
                  h(layer) = q(i,j,m+1)
                  if (h(layer) > drytolerance) then
                      u(layer) = q(i,j,m+2) / q(i,j,m+1)
                      v(layer) = q(i,j,m+3)/ q(i,j,m+1)
                  else
                      dry_state(layer) = .true.
                      u(layer) = 0.d0
                      v(layer) = 0.d0
                  endif
              enddo
              if (sum(h) > drytolerance) then
                  aux(i,j,9)=(u(1) - u(2))**2 / (g*one_minus_r*sum(h))
                  if ((aux(i,j,9) > richardson_tolerance)
     &                  .and.(.not.dry_state(2))) then
                      print 100,i,j,aux(i,j,9)
                  endif
                  aux(i,j,10)=(v(1) - v(2))**2 / (g*one_minus_r*sum(h))
                  if ((aux(i,j,10) > richardson_tolerance)
     &                  .and.(.not.dry_state(2))) then
                      print 100,i,j,aux(i,j,10)
                  endif
               else
                   aux(i,j,9) = 0.d0
                   aux(i,j,10) = 0.d0
               endif
          enddo
      enddo
      endif
      
100   format ("Hyperbolicity failed, kappa(",i4,",",i4,") = ",d16.8)

    
C     ! These need to be modified to use other aux array locations, 7,8 are 
C     ! being used for initial etas
      ! Calculate gradient of Pressure
C       aux(:,:,7) = 0.d0
C       aux(:,:,8) = 0.d0
C       do i=1,mx
C          do j=1,my
C                if (abs(q(i,j,1)) > drytolerance) then
C                    diff_x = aux(i+1,j,6) - aux(i-1,j,6)
C                    diff_y = aux(i,j+1,6) - aux(i,j-1,6)
C                    if (abs(diff_x) < pressure_tolerance) then
C                        aux(i,j,7) = 0.d0
C                    else
C                        aux(i,j,7) = (diff_x) / (2.d0 * 1000.d0 * dx)
C                    endif
C                    if (abs(diff_y) < pressure_tolerance) then
C                        aux(i,j,8) = 0.d0
C                    else
C                        aux(i,j,8) = (diff_y) / (2.d0 * 1000.d0 * dy)
C                    endif   
C                endif
C          enddo
C       enddo
C       
C       ! Calculate vorticity
C       aux(:,:,9) = 0.d0
C       vel = 0.d0
C       do i=1,mx
C           do j=1,my
C               if (abs(q(i,j,1)) > drytolerance) then
C                   vel(i,j,1) = q(i,j,2) / q(i,j,1)
C                   vel(i,j,2) = q(i,j,3) / q(i,j,1)
C               endif
C           enddo
C       enddo
C       do i=1,mx
C           do j=1,my
C               aux(i,j,9) = (vel(i+1,j,2) - vel(i-1,j,2)) / (2.d0*dx) -         
C      &                     (vel(i,j+1,1) - vel(i,j-1,1)) / (2.d0*dy)
C           enddo
C       enddo
      
      return
      end
    