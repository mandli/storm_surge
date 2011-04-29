subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

    use hurricane_module
    use geoclaw_module

    implicit none
    
    ! Input parameters
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Ouput
    double precision, intent(inout) :: q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
    double precision, intent(inout) :: aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
    
    ! Locals
    integer :: i,j
    double precision :: h,g,coeff,hu,hv,gamma,dgamma,tol,w,ycell,cor,hu0,hv0
    double precision :: a11,a12,a21,a22,ct,xc,yc,cd,tau,wind_speed
    double precision :: P_atmos_x,P_atmos_y
    
    double precision, parameter :: rho_water = 1d3
    
    integer :: mn,n
    double precision :: local_dt,local_t

    ! Check for NANs in solution:
!     call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,2)

    g=grav
    coeff = coeffmanning
    tol = 1.d-30  ! To prevent divide by zero in gamma

    ! friction--------------------------------------------------------
    if (coeffmanning.gt.0.d0.and.frictiondepth.gt.0.d0) then
        do i=1,mx
            do j=1,my
                h=q(i,j,1)
                if (h.le.frictiondepth) then
                    ! Apply friction source term only in shallower water
                    hu=q(i,j,2)
                    hv=q(i,j,3)
    
                    if (h.lt.tol) then
                        q(i,j,2)=0.d0
                        q(i,j,3)=0.d0
                    else
                        gamma= dsqrt(hu**2 + hv**2)*(g*coeff**2)/(h**(7/3))
                        dgamma=1.d0 + dt*gamma
                        q(i,j,2)= q(i,j,2)/dgamma
                        q(i,j,3)= q(i,j,3)/dgamma
                    endif
                endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! coriolis--------------------------------------------------------
    if (icoordsys.eq.2.and.icoriolis.eq.1) then
        w = 2.d0*pi/(86400.d0) !angular velocity of earth
        do i=1,mx
            do j=1,my
                ycell = ylower + (j-.5d0)*dy
                cor = 2.d0*w*sin(pi*ycell/180.d0)
                ct = cor*dt
                !integrate momentum exactly using matrix exponential
                !forth order term should be sufficient since cor^3 ~= eps
                hu0 = q(i,j,2)
                hv0 = q(i,j,3)
                !dq/dt = 2w*sin(latitude)*[0 1 ; -1 0] q = Aq
                !e^Adt = [a11 a12; a21 a22] + I
                a11 = -0.5d0*ct**2 + ct**4/24.d0
                a12 = ct - ct**3/6.0d0
                a21 = -ct + ct**3/6.0d0
                a22 = a11
                !q = e^Adt * q0
                q(i,j,2) = q(i,j,2) + hu0*a11 + hv0*a12
                q(i,j,3) = q(i,j,3) + hu0*a21 + hv0*a22
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        do i=1,mx
            do j=1,my
!                 if (abs(q(i,j,1)) > drytolerance) then
                    wind_speed = sqrt(aux(i,j,4)**2 + aux(i,j,5)**2)
                    tau = wind_drag(wind_speed) * rho_air * wind_speed
                    q(i,j,2) = q(i,j,2) + dt * tau * aux(i,j,4) / rho_water
                    q(i,j,3) = q(i,j,3) + dt * tau * aux(i,j,5) / rho_water
!                 endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! atmosphere -----------------------------------------------------
    if (pressure_forcing) then
        do i=1,mx
            do j=1,my                    
                q(i,j,2) = q(i,j,2) - dt * aux(i,j,7) * q(i,j,1)
                q(i,j,3) = q(i,j,3) - dt * aux(i,j,8) * q(i,j,1)
            enddo
        enddo
                
    endif
    ! ----------------------------------------------------------------
end subroutine src2