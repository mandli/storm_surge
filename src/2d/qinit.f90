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
    double precision, intent(inout) :: q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc)
    
    ! Locals
    integer :: i,j
!     double precision :: x,y,xim,xip,yjm,yjp,xc,yc,dq,total_depth
    double precision :: x,y,xmid,g,m,x_c,y_c
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
            
            q(1,i,j) = aux(7,i,j) * rho(1)
            q(2:3,i,j) = 0.d0
            if (layers > 1) then
                q(4,i,j) = aux(8,i,j) * rho(2)
                q(5:6,i,j) = 0.d0
            endif
            
            ! Test perturbations - these only work in the x-direction
            if (init_type == 1 .or. init_type == 2) then
                ! Calculate wave family for perturbation
                gamma = aux(8,i,j) / aux(7,i,j)
                select case(wave_family)
                    case(1) ! Shallow water, left-going
                        alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                        lambda = -sqrt(g*aux(7,i,j)*(1.d0+alpha))
                    case(2) ! Internal wave, left-going
                        alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                        lambda = -sqrt(g*aux(7,i,j)*(1.d0+alpha))
                    case(3) ! Internal wave, right-going
                        alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                        lambda = sqrt(g*aux(7,i,j)*(1.d0+alpha))
                    case(4) ! Shallow water, right-going
                        alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                        lambda = sqrt(g*aux(7,i,j)*(1.d0+alpha))
                end select
                eigen_vector = [1.d0,lambda,0.d0,alpha,lambda*alpha,0.d0]

                if (init_type == 1) then
                    ! Add perturbation
                    if ((x < init_location(1)).and.(wave_family >= 3)) then
                        q(1:3,i,j) = q(1:3,i,j) + rho(1) * epsilon * eigen_vector(1:3)
                        q(4:6,i,j) = q(4:6,i,j) + rho(2) * epsilon * eigen_vector(4:6)
                    else if ((x > init_location(1)).and.(wave_family < 3)) then
                        q(1:2,i,j) = q(1:2,i,j) + rho(1) * epsilon * eigen_vector(1:2)
                        q(4:5,i,j) = q(4:5,i,j) + rho(2) * epsilon * eigen_vector(4:5)
                    endif
                ! Gaussian wave along a direction on requested wave family
                else if (init_type == 2) then
                    ! Transform back to computational coordinates
                    x_c = x * cos(angle) + y * sin(angle) - init_location(1)
                    deta = epsilon * exp(-(x_c/sigma)**2)
                    q(1,i,j) = q(1,i,j) + rho(1) * deta
                    q(4,i,j) = q(4,i,j) + rho(2) * alpha * deta
                endif
            ! Symmetric gaussian hump
            else if (init_type == 3) then
                deta = epsilon * exp(-((x-init_location(1))/sigma)**2)  &
                               * exp(-((y-init_location(2))/sigma)**2)
                q(1,i,j) = q(1,i,j) + rho(1) * deta
            ! Shelf conditions from AN paper
            else if (init_type == 4) then
                alpha = 0.d0
                xmid = 0.5d0*(-180.e3 - 80.e3)
                if ((x > -130.e3).and.(x < -80.e3)) then
                    deta = epsilon * sin((x-xmid)*PI/(-80.e3-xmid))
                    q(4,i,j) = q(4,i,j) + rho(2) * alpha * deta
                    q(1,i,j) = q(1,i,j) + rho(1) * deta * (1.d0 - alpha)
                endif
            ! Inundation test
            else if (init_type == 5) then
                x_c = (x - init_location(1)) * cos(angle) &
                    + (y - init_location(2)) * sin(angle)
                deta = epsilon * exp(-(x_c/sigma)**2)
                q(1,i,j) = q(1,i,j) + rho(1) * deta
            endif
            
        enddo
    enddo

end subroutine qinit
