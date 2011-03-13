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
    double precision :: x,eigen_vector(4),gamma,lambda,alpha
    
    do i=1,mx
        x = xlower+(i-0.5)*dx

        ! Set initial perturbation to zero
        q(i,:) = 0.d0
        
        ! Calculate wave family for perturbation
        gamma = aux(i,3) / aux(i,2)
        select case(wave_family)
            case(1) ! Shallow water, left-going
                alpha = 0.5d0 * (gamma - 1 + sqrt((gamma-1)**2+4.d0*r*gamma))
                lambda = -sqrt(g*aux(i,2)*(1+alpha))
            case(2) ! Internal wave, left-going
                alpha = 0.5d0 * (gamma - 1 - sqrt((gamma-1)**2+4.d0*r*gamma))
                lambda = -sqrt(g*aux(i,2)*(1+alpha))
            case(3) ! Internal wave, right-going
                alpha = 0.5d0 * (gamma - 1 - sqrt((gamma-1)**2+4.d0*r*gamma))
                lambda = sqrt(g*aux(i,2)*(1+alpha))
            case(4) ! Shallow water, right-going
                alpha = 0.5d0 * (gamma - 1 + sqrt((gamma-1)**2+4.d0*r*gamma))
                lambda = sqrt(g*aux(i,2)*(1+alpha))
        end select
        eigen_vector = [1.d0,lambda,alpha,lambda*alpha]
            
        ! Add perturbation
        if (x < init_location) then
            q(i,:) = q(i,:) + epsilon * eigen_vector
        endif
    enddo
    
end subroutine qinit