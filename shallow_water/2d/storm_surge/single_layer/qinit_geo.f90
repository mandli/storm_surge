subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use multilayer_module
    use geoclaw_module
    use hurricane_module

    implicit none
    
    ! Input
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy
    double precision, intent(in) :: aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

    ! Input/Output
    double precision, intent(inout) :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

    ! Locals
    integer :: i,j
    double precision :: x,y,wind_speed,tau,deta

    do i=1-mbc,mx+mbc
        x = xlower + (i-0.5d0)*dx
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0)*dy
            q(i,j,1)=dmax1(0.d0,sealevel-aux(i,j,1))
            q(i,j,2)=0.d0
            q(i,j,3)=0.d0
            
            if (init_type == 3) then
                deta = epsilon * exp(-((x-init_location(1))/sigma)**2) * exp(-((y-init_location(2))/sigma)**2)
                q(i,j,1) = q(i,j,1) + deta
            else if (init_type == 7) then
                if (aux(i+1,j,1) > -300.d0 .and. .not.(aux(i,j,1) > -300.d0)) then
                    q(i,j,1) = q(i,j,1) + epsilon
                endif
            endif
        enddo
    enddo
    
    ! Get and store hurricane wind field
    call hurricane_wind(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,-ramp_up_time,q(:,:,2:3))
      
    ! Setup swirling velocity field
!     do i=1-mbc,mx+mbc
!         do j=1-mbc,my+mbc
! 
!             ! Scale with appropriate factors
!             wind_speed = sqrt(q(i,j,2)**2 + q(i,j,3)**2)
!             tau = wind_drag(wind_speed) * rho_air * wind_speed
!             q(i,j,2) = tau * q(i,j,2) * q(i,j,1)
!             q(i,j,3) = tau * q(i,j,3) * q(i,j,1)
!         enddo
!     enddo
end subroutine qinit