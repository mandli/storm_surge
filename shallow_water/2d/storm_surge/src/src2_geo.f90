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
    double precision :: yc
    ! Friction
    double precision :: h(2),g,coeff,tol,speed,D
    ! Wind and pressure
    double precision :: tau,wind_speed,P_atmos_x,P_atmos_y
    ! Coriolis source term
    ! angular velocity of earth = 2.d0*pi/(86400.d0) 
!     double precision, parameter :: OMEGA = 7.2722052166430395d-05
!     double precision, parameter :: OMEGA = 7.2722052166430395d-03 ! Incorrect value but used for testing
    ! For beta-plane approximation, center of domain's latitude 
    double precision :: fdt,theta,a11,a12,a21,a22,hu,hv
    
    integer :: mn,n

    ! Check for NANs in solution:
!     call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,2)

    g=grav
    coeff = coeffmanning
    tol = 1.d-30  ! To prevent divide by zero in gamma

    ! friction--------------------------------------------------------
    if (coeffmanning > 0.d0 .and. ifriction > 1) then
        ! Constant coefficient of friction
        if (ifriction == 1) then
            D = coeff
            do i=1,mx
                do j=1,my
                    ! Check to see which layer we are doing this on
                    if (layers > 1) then
                        h(1) = q(i,j,1) / rho(1)
                        h(2) = q(i,j,4) / rho(2)
                    else
                        h(1) = q(i,j,1)
                        h(2) = 0.d0
                    endif
            
                    ! Bottom layer wet, apply to bottom layer
                    if (h(2) > tol) then
                        ! Exactly integrate and modify bottom layer momentum
                        q(i,j,5) = q(i,j,5) * exp(-D*dt)
                        q(i,j,6) = q(i,j,6) * exp(-D*dt)
                
                    ! Only top layer wet, apply to top layer only
                    else if (h(1) > tol) then
                        ! Set bottom layer momentum to zero
                        if (layers > 1) q(i,j,5:6) = 0.d0
                        
                        ! Exactly integrate and modify top layer momentum
                        q(i,j,2) = q(i,j,2) * exp(-D*dt)
                        q(i,j,3) = q(i,j,3) * exp(-D*dt)
                
                    ! Neither layer wet, set momentum to zero
                    else
                        q(i,j,2:3) = 0.d0
                        if (layers > 1) q(i,j,5:6) = 0.d0
                    endif
                enddo
            enddo
        ! Manning-N friction
        else if (ifriction == 2) then
            do i=1,mx
                do j=1,my
                    ! Check to see which layer we are doing this on
                    if (layers > 1) then
                        h(1) = q(i,j,1) / rho(1)
                        h(2) = q(i,j,4) / rho(2)
                    else
                        h(1) = q(i,j,1)
                        h(2) = 0.d0
                    endif
            
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
                        if (layers > 1) q(i,j,5:6) = 0.d0
                
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
                        if (layers > 1) q(i,j,5:6) = 0.d0
                    endif
                enddo
            enddo
        endif
    endif
    ! ----------------------------------------------------------------

    ! coriolis--------------------------------------------------------
    if (icoriolis > 0) then
        do j=1,my
            yc = ylower + (j-.5d0)*dy
            fdt = coriolis(yc) * dt
            do i=1,mx
                !dq/dt = 2w*sin(latitude)*[0 1 ; -1 0] q = Aq
                !e^Adt = [a11 a12; a21 a22] + I
                a11 = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
                a12 = fdt - fdt**3 / 6.0d0
                a21 = -fdt + fdt**3 / 6.0d0
                a22 = a11
                hu = q(i,j,2)
                hv = q(i,j,3)
                !q = e^Adt * q0
                q(i,j,2) = hu * a11 + hv * a12
                q(i,j,3) = hu * a21 + hv * a22
                if (layers > 1) then
                    hu = q(i,j,5)
                    hv = q(i,j,6)
                    q(i,j,5) = hu * a11 + hv * a12
                    q(i,j,6) = hu * a21 + hv * a22
                endif
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
                h = 0.d0  
                if (layers > 1) then
                    h(1) = q(i,j,1) / rho(1)
                    h(2) = q(i,j,4) / rho(2)
                else
                    h(1) = q(i,j,1)
                endif
                
                ! Calculate gradient of Pressure
                P_atmos_x = (aux(i+1,j,6) - aux(i-1,j,6)) / (2.d0*dx)
                P_atmos_y = (aux(i,j+1,6) - aux(i,j-1,6)) / (2.d0*dy)
                if (abs(P_atmos_x) < pressure_tolerance) then
                    P_atmos_x = 0.d0
                endif
                if (abs(P_atmos_y) < pressure_tolerance) then
                    P_atmos_y = 0.d0
                endif
                
                if (layers > 1) then
                    if (h(1) > drytolerance) then
                        q(i,j,2) = q(i,j,2) - dt * h(1) * P_atmos_x
                        q(i,j,3) = q(i,j,3) - dt * h(1) * P_atmos_y
                    else if (h(2) > drytolerance) then
                        q(i,j,5) = q(i,j,5) - dt * h(2) * P_atmos_x
                        q(i,j,6) = q(i,j,6) - dt * h(2) * P_atmos_y
                    endif
                else
                    if (h(1) > drytolerance) then
                        q(i,j,2) = q(i,j,2) - dt * h(1) * P_atmos_x / rho(1)
                        q(i,j,3) = q(i,j,3) - dt * h(1) * P_atmos_y / rho(1)
                    endif
                endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------
end subroutine src2