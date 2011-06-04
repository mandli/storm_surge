! This routine should be a simplified version of src2
! which applies source terms for a 1-d slice of data along the
! edge of a grid.  This is called only from qad where the conservative
! fix-up is applied and is used to apply source terms over partial
! time steps to the coarse grid cell values used in solving Riemann
! problems at the interface between coarse and fine grids.
subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)

    use multilayer_module, only: rho,eta,layers
    use hurricane_module
    use geoclaw_module

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx1d,maux
    double precision, intent(in) :: t,dt
    
    ! Output
    double precision, intent(inout) :: q1d(mx1d,meqn)
    double precision, intent(inout) :: aux1d(mx1d,maux)

    ! Locals
    integer :: i
    double precision :: g,coeff,tol,h(2)
    double precision :: wind_speed,tau,P_atmos_x,P_atmos_y,speed,D

    ! Common block
    double precision dtcom,dxcom,dycom,tcom
    integer icom,jcom
    
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

    ! Incorporates friction using Manning coefficient
    g=grav
    coeff = coeffmanning
    tol = 1.d-30  !# to prevent divide by zero in gamma

    ! = Friction =============================================================
    if (coeff > 0.d0) then
        do i=1,mx1d
            h(1) = q1d(i,1) / rho(1)
            if (layers > 1) then
                h(2) = q1d(i,4) / rho(2)
            else
                h(2) = 0.d0
            endif
            
            if (h(2) > tol) then
                speed = sqrt(q1d(i,5)**2 + q1d(i,6)**2) / q1d(i,4)
                D = coeff**2 * g * sum(h)**(-7/3) * speed
                
                q1d(i,5) = q1d(i,5) * exp(-D*dt)
                q1d(i,6) = q1d(i,6) * exp(-D*dt)
            else if (h(1) > tol) then
                if (layers > 1) q1d(i,5:6) = 0.d0
                
                speed = sqrt(q1d(i,2)**2 + q1d(i,3)**2) / q1d(i,1)
                
                D = coeff**2 * g * sum(h)**(-7/3) * speed
                
                q1d(i,2) = q1d(i,2) * exp(-D*dt)
                q1d(i,3) = q1d(i,3) * exp(-D*dt)
                
            else
                q1d(i,2:3) = 0.d0
                if (layers > 1)  q1d(i,5:6) = 0.d0
            endif
        enddo          
    endif
    
    ! = Wind Forcing =========================================================
    if (wind_forcing) then
        ! Here we have to take into account geoclaw which we use for the 
        ! single layer case.  It needs to divide by the water's density
        if (layers > 1) then
            do i=1,mx1d
                if (abs(q1d(i,1)) > drytolerance) then
                    wind_speed = sqrt(aux1d(i,4)**2 + aux1d(i,5)**2)
                    if (wind_speed > wind_tolerance) then
                        tau = wind_drag(wind_speed) * rho_air * wind_speed
                        q1d(i,2) = q1d(i,2) + dt * tau * aux1d(i,4)
                        q1d(i,3) = q1d(i,3) + dt * tau * aux1d(i,5)
                    endif
                endif
            enddo
        else
            do i=1,mx1d
                if (abs(q1d(i,1)) > drytolerance) then
                    wind_speed = sqrt(aux1d(i,4)**2 + aux1d(i,5)**2)
                    tau = wind_drag(wind_speed) * rho_air * wind_speed
                    q1d(i,2) = q1d(i,2) + dt * tau * aux1d(i,4) / rho(1)
                    q1d(i,3) = q1d(i,3) + dt * tau * aux1d(i,5) / rho(1)
                endif
            enddo
        endif
    endif
    ! ========================================================================
    
    ! == Pressure Forcing ====================================================
    if (pressure_forcing) then
        stop "Not sure how to proceed, need direction or the right dx or dy."
    endif

end subroutine src1d
