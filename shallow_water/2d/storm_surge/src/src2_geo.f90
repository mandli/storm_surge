subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

    use multilayer_module, only: rho,eta,layers
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
    double precision :: w,ycell,cor,hu0,hv0
    double precision :: a11,a12,a21,a22,ct,xc,yc,cd
    double precision :: h(2),g,coeff,tol,speed,D
    double precision :: tau,wind_speed,P_atmos_x,P_atmos_y
    
    integer :: mn,n

    ! Check for NANs in solution:
!     call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,2)

    g=grav
    coeff = coeffmanning
    tol = 1.d-30  ! To prevent divide by zero in gamma

    ! friction--------------------------------------------------------
    if (coeffmanning > 0.d0) then
        do i=1,mx
            do j=1,my
                ! Check to see which layer we are doing this on
                h(1) = q(i,j,1) / rho(1)
                h(2) = q(i,j,4) / rho(2)
                
                ! Bottom layer wet, apply to bottom layer
                if (h(2) > tol) then
                    ! Extract speed of bottom layer
                    speed = sqrt(q(i,j,5)**2 + q(i,j,6)**2) / q(i,j,4)
                                        
                    ! Calculate drag coefficient 
                    D = coeff**2 * g * sum(h)**(-7/3) * speed
                    
                    ! Exactly integrate and modify bottom layer momentum
                    q(i,j,5) = q(i,j,5) * exp(-D*dt)
                    q(i,j,6) = q(i,j,6) * exp(-D*dt)
                    
                ! Only top layer wet, apply to top layer only
                else if (h(1) > tol) then
                    ! Set bottom layer momentum to zero
                    q(i,j,5:6) = 0.d0
                    
                    ! Extract speed of top layer
                    speed = sqrt(q(i,j,2)**2 + q(i,j,3)**2) / q(i,j,1)
                    
                    ! Calculate drag coefficient
                    D = coeff**2 * g * sum(h)**(-7/3) * speed
                    
                    ! Exactly integrate and modify top layer momentum
                    q(i,j,2) = q(i,j,2) * exp(-D*dt)
                    q(i,j,3) = q(i,j,3) * exp(-D*dt)
                    
                ! Neither layer wet, set momentum to zero
                else
                    q(i,j,2:3) = 0.d0
                    q(i,j,5:6) = 0.d0
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
        ! Here we have to take into account geoclaw which we use for the 
        ! single layer case.  It needs to divide by the water's density
        if (layers > 1) then
            do i=1,mx
                do j=1,my
                    if (q(i,j,1) / rho(1) > drytolerance) then
                        wind_speed = sqrt(aux(i,j,4)**2 + aux(i,j,5)**2)
                        if (wind_speed > wind_tolerance) then
                            tau = wind_drag(wind_speed) * rho_air * wind_speed
                            q(i,j,2) = q(i,j,2) + dt * tau * aux(i,j,4)
                            q(i,j,3) = q(i,j,3) + dt * tau * aux(i,j,5)
                        endif
                    endif
                enddo
            enddo
        else
            do i=1,mx
                do j=1,my
    !                 if (abs(q(i,j,1)) > drytolerance) then
                        wind_speed = sqrt(aux(i,j,4)**2 + aux(i,j,5)**2)
                        tau = wind_drag(wind_speed) * rho_air * wind_speed
                        q(i,j,2) = q(i,j,2) + dt * tau * aux(i,j,4) / rho(1)
                        q(i,j,3) = q(i,j,3) + dt * tau * aux(i,j,5) / rho(1)
    !                 endif
                enddo
            enddo
        endif
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