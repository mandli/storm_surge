! This routine should be a simplified version of src2
! which applies source terms for a 1-d slice of data along the
! edge of a grid.  This is called only from qad where the conservative
! fix-up is applied and is used to apply source terms over partial
! time steps to the coarse grid cell values used in solving Riemann
! problems at the interface between coarse and fine grids.
subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)

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
    double precision :: g,coeff,tol,h,hu,hv,gamma,dgamma
    double precision :: wind_speed,tau,P_atmos_x,P_atmos_y
    
    double precision, parameter :: rho_water = 1d3

    ! Common block
    double precision dtcom,dxcom,dycom,tcom
    integer icom,jcom
    
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

    ! Incorporates friction using Manning coefficient
    g=grav
    coeff = coeffmanning
    tol = 1.d-30  !# to prevent divide by zero in gamma

    if ((coeffmanning > 0.d0).and.(frictiondepth > 0.d0)) then
        do i=1,mx1d
            h=q1d(i,1)
            if (h.lt.frictiondepth) then
            ! Apply friction source term only in shallower water
                hu=q1d(i,2)
                hv=q1d(i,3)

                if (h.lt.tol) then
                    q1d(i,2)=0.d0
                    q1d(i,3)=0.d0
                else
                    gamma= dsqrt(hu**2 + hv**2)*(g*coeff**2)/(h**(7/3))
                    dgamma=1.d0 + dt*gamma
                    q1d(i,2)= q1d(i,2)/dgamma
                    q1d(i,3)= q1d(i,3)/dgamma
                endif
            endif
        enddo
    endif
    
    if (wind_forcing) then
        do i=1,mx1d
            if (abs(q1d(i,1)) > drytolerance) then
                wind_speed = sqrt(aux1d(i,4)**2 + aux1d(i,5)**2)
                tau = wind_drag(wind_speed) * rho_air * wind_speed
                q1d(i,2) = q1d(i,2) + dt * tau * aux1d(i,4) / rho_water
                q1d(i,3) = q1d(i,3) + dt * tau * aux1d(i,5) / rho_water
            endif
        enddo
    endif
    
    if (pressure_forcing) then
        do i=1,mx1d
            q1d(i,2) = q1d(i,2) - dt * aux1d(i,7) * q1d(i,1)
            q1d(i,3) = q1d(i,3) - dt * aux1d(i,8) * q1d(i,1)
        enddo
    endif

end subroutine src1d
