subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)

    use parameters_module

    implicit none
    
    ! Input arguments
    integer, intent(in) :: maxmx,meqn,mbc,mx,maux
    double precision, intent(in) :: xlower,dx
    
    double precision, intent(inout) :: q(1-mbc:maxmx+mbc,meqn)
    double precision, intent(inout) :: aux(1-mbc:maxmx+mbc,maux)

    ! Locals
    integer :: i
    double precision :: x,eigen_vector(4),gamma,lambda,alpha,h_1,h_2

    if (.not.((0 <= init_type).and.(init_type <= 3))) then
        print "(a,i2)","Invalid initialization type requested, init_type = ",init_type
    endif
    
    do i=1,mx
        x = xlower+(i-0.5)*dx

        ! Set initial perturbation to zero
        q(i,1) = aux(i,3) * rho(1)
        q(i,3) = aux(i,4) * rho(2)
        q(i,2) = 0.d0
        q(i,4) = 0.d0
        
        if (init_type == 1) then
            ! Calculate wave family for perturbation
            gamma = aux(i,4) / aux(i,3)
            select case(wave_family)
                case(1) ! Shallow water, left-going
                    alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                    lambda = -sqrt(g*aux(i,3)*(1.d0+alpha))
                case(2) ! Internal wave, left-going
                    alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                    lambda = -sqrt(g*aux(i,3)*(1.d0+alpha))
                case(3) ! Internal wave, right-going
                    alpha = 0.5d0 * (gamma - 1.d0 - sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                    lambda = sqrt(g*aux(i,3)*(1.d0+alpha))
                case(4) ! Shallow water, right-going
                    alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
                    lambda = sqrt(g*aux(i,3)*(1.d0+alpha))
            end select
            eigen_vector = [1.d0,lambda,alpha,lambda*alpha]
            
            ! Add perturbation
            if ((x < init_location).and.(wave_family >= 3)) then
                q(i,1:2) = q(i,1:2) + rho(1) * epsilon * eigen_vector(1:2)
                q(i,3:4) = q(i,3:4) + rho(2) * epsilon * eigen_vector(3:4)
            else if ((x > init_location).and.(wave_family < 3)) then
                q(i,1:2) = q(i,1:2) + rho(1) * epsilon * eigen_vector(1:2)
                q(i,3:4) = q(i,3:4) + rho(2) * epsilon * eigen_vector(3:4)
            endif
        else if (init_type == 2) then
            gamma = aux(i,4) / aux(i,3)
            alpha = 0.5d0 * (gamma - 1.d0 + sqrt((gamma-1.d0)**2+4.d0*r*gamma))
            q(i,1) = q(i,1) + rho(1) * epsilon * exp(-((x-init_location)/sigma)**2)
            q(i,3) = q(i,3) + rho(2) * alpha * epsilon * exp(-((x-init_location)/sigma)**2)
        else if (init_type == 3) then
            q(i,1) = q(i,1) - rho(1) * epsilon * exp(-((x-init_location)/sigma)**2)
            q(i,3) = q(i,3) + rho(2) * epsilon * exp(-((x-init_location)/sigma)**2)
        endif
    enddo
    
end subroutine qinit