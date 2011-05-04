subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
!     =====================================================
!
!      # Set initial sea level flat unless iqinit = 1, in which case
!      # an initial perturbation of the q(i,j,1) is specified and has
!      # been strored in qinitwork.

    use multilayer_module
    use geoclaw_module

    implicit none

    ! Input arguments
    integer, intent(in) :: maxmx,maxmy,meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy
    double precision, intent(inout) :: q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
    double precision, intent(inout) :: aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
    
    ! Locals
    integer :: i,j
!     double precision :: x,y,xim,xip,yjm,yjp,xc,yc,dq,total_depth
    double precision :: x,y,xmid,g,m
    double precision :: eigen_vector(6),gamma,lambda,alpha,h_1,h_2,deta

    g = grav
    ! External definitions
!     include 'qinit.i'
!     external topointegral

    ! Set the initial state to be flat
    do i=1,mx
        x = xlower + (i-0.5d0)*dx
        do j=1,my
            y = ylower + (j-0.5d0)*dy
            
            q(i,j,1) = aux(i,j,7) * rho(1)
            q(i,j,2:3) = 0.d0
            if (layers > 1) then
                q(i,j,4) = aux(i,j,8) * rho(2)
                q(i,j,5:6) = 0.d0
            endif
            
            ! Test perturbations - these only work in the x-direction
            if (init_type == 1) then
                ! Calculate wave family for perturbation
                gamma = aux(i,j,8) / aux(i,j,7)
                select case(wave_family)
                    case(1) ! Shallow water, left-going
                        alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                        lambda = -sqrt(g*aux(i,j,7)*(1.d0+alpha))
                    case(2) ! Internal wave, left-going
                        alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                        lambda = -sqrt(g*aux(i,j,7)*(1.d0+alpha))
                    case(3) ! Internal wave, right-going
                        alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                        lambda = sqrt(g*aux(i,j,7)*(1.d0+alpha))
                    case(4) ! Shallow water, right-going
                        alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                        lambda = sqrt(g*aux(i,j,7)*(1.d0+alpha))
                end select
                eigen_vector = [1.d0,lambda,0.d0,alpha,lambda*alpha,0.d0]
            
                ! Add perturbation
                if ((x < init_location(1)).and.(wave_family >= 3)) then
                    q(i,j,1:3) = q(i,j,1:3) + rho(1) * epsilon * eigen_vector(1:3)
                    q(i,j,4:6) = q(i,j,4:6) + rho(2) * epsilon * eigen_vector(4:6)
                else if ((x > init_location(1)).and.(wave_family < 3)) then
                    q(i,j,1:2) = q(i,j,1:2) + rho(1) * epsilon * eigen_vector(1:2)
                    q(i,j,4:5) = q(i,j,4:5) + rho(2) * epsilon * eigen_vector(4:5)
                endif
            else if (init_type == 2) then
                gamma = aux(i,j,8) / aux(i,j,7)
                alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                deta = epsilon * exp(-((x-init_location(1))/sigma)**2)
                q(i,j,1) = q(i,j,1) + rho(1) * deta
                q(i,j,4) = q(i,j,4) + rho(2) * alpha * deta
            ! Gaussian hump on top surface
            elseif (init_type == 3) then
                deta = epsilon * exp(-((x-init_location(1))/sigma)**2) * exp(-((y-init_location(2))/sigma)**2)
                q(i,j,1) = q(i,j,1) + rho(1) * deta
            ! Gassian hump on internal surface
            else if (init_type == 4) then
                deta = epsilon * exp(-((x-init_location(1))/sigma)**2) * exp(-((y-init_location(2))/sigma)**2)
                q(i,j,1) = q(i,j,1) - rho(1) * deta
                q(i,j,4) = q(i,j,4) + rho(2) * deta
            else if (init_type == 5) then
                gamma = aux(i,j,8) / aux(i,j,7)
    !             alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                alpha = 0.d0
                xmid = 0.5d0*(-180.e3+-80.e3)
!                 xmid = 300d3
                if ((x > 275.d3).and.(x < 325.d3)) then
                    deta = epsilon * sin((x-xmid)*PI/(-80.e3-xmid))
                    q(i,j,4) = q(i,j,4) + rho(2) * alpha * deta
                    q(i,j,1) = q(i,j,1) + rho(1) * deta * (1.d0 - alpha)
                endif
            else if (init_type == 6) then
                gamma = aux(i,j,8) / aux(i,j,7)
                alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                deta = epsilon * exp(-((y-init_location(2))/sigma)**2)
                q(i,j,1) = q(i,j,1) + rho(1) * deta
                q(i,j,4) = q(i,j,4) + rho(2) * alpha * deta
            else if (init_type == 7) then
                if (eta(2) < aux(i+1,j,1) .and. .not.(eta(2) < aux(i,j,1))) then
!                 if (q(i+1,j,4) / rho(2) < drytolerance .and. .not.(q(i,j,4) / rho(2) < drytolerance)) then
                    q(i,j,1) = q(i,j,1) + epsilon * rho(1)
!                     q(i,j,4) = q(i,j,4) + epsilon * rho(2)
                endif
            else if (init_type == 8) then
                if (init_location(1) < x) then
                    aux(i,j,7) = eta(1) - aux(i,j,1)
                    aux(i,j,8) = 0.d0
                    q(i,j,1) = aux(i,j,7) * rho(1)
                    q(i,j,4) = aux(i,j,8) * rho(2)
                endif
            endif
        enddo
    enddo

end subroutine qinit
