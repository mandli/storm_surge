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
    double precision, dimension(2) :: u_l,h_l,u_r,h_r
    double precision :: x_center,eig_vec(4),lambda,mult_depth,location

    ! Set base state
    q(:,1) = depth_ratio * total_depth
    q(:,2) = 0.d0
    q(:,3) = total_depth - q(:,1)
    q(:,4) = 0.d0
    
    select case(init_type)
        ! Constant initial conditions
        case(0)
            
        ! Riemann initial conditions
        case(1)
            ! Left state is assumed to be base state
            forall(i=1:mx,xlower+(i-0.5d0)*dx < init_location)
                q(i,1) = 0.75d0           ! h1
                q(i,2) = 0.d0            ! h1 u1
                q(i,3) = 0.3d0           ! h2
                q(i,4) = 0.d0            ! h2 u2
            end forall
    
            ! Right state
            forall(i=1:mx,xlower+(i-0.5d0)*dx >= init_location)
                q(i,1) = 0.75d0           ! h1
                q(i,2) = 0.d0            ! h1 u1
                q(i,3) = 0.25d0           ! h2
                q(i,4) = 0.d0            ! h2 u2
            end forall
            
        ! Smooth gaussian hump of water located at jump_location
        case(2)
            ! Base state
            q(:,1) = 0.75d0
            q(:,3) = 0.25d0

            forall(i=1:mx)
                q(i,3) = q(i,3) + beta * exp(-((xlower+(i-0.5d0)*dx-init_location)/sigma)**2)
            end forall
        
        ! Perturbation in one of the wave families    
        case(3)
            ! Add perturbation of eigenvector to initial condition
            ! Currently for 1st and 4th wave families
            eig_vec(1) = 1.d0
            one_minus_r = 1.d0 - r
            do i=1,mx
                if (xlower+(i-0.5d0)*dx < init_location) then
                    mult_depth = q(i,1) * q(i,3)
                    lambda = - sqrt(g*total_depth) + 0.5d0 * mult_depth/ total_depth**(3/2) * one_minus_r
!                     lambda = sqrt(g*total_depth) + 0.5d0 * mult_depth/ total_depth**(3/2) * one_minus_r
                    eig_vec(2) = lambda
                    eig_vec(3) = (lambda**2 - g*q(i,1)) / (g*r*q(i,1))
                    eig_vec(4) = eig_vec(3) * eig_vec(2)
                    q(i,:) = q(i,:) + beta * eig_vec
                endif
            enddo
        case(4)
            ! Customized init
            h_l = [0.5d0,0.5d0]
            u_l = 0.d0
            h_r = [0.0d0,0.3d0]
            u_r = 0.d0
            
            ! Dry state Riemann problem
            forall(i=1:mx,xlower+(i-0.5d0)*dx >= init_location)
                q(i,1) = h_r(1)           ! h1
                q(i,2) = h_r(1)*u_r(1)    ! h1 u1
                q(i,3) = h_r(2)           ! h2
                q(i,4) = h_r(2)*u_r(2)    ! h2 u2
            end forall
            forall(i=1:mx,xlower+(i-0.5d0)*dx < init_location)
                q(i,1) = h_l(1)           ! h1
                q(i,2) = h_l(1)*u_l(1)    ! h1 u1
                q(i,3) = h_l(2)           ! h2
                q(i,4) = h_l(2)*u_l(2)    ! h2 u2
            end forall

            ! Smooth gaussian hump
            sigma = 0.02d0
            beta = 0.05d0
            location = -0.25d0
            forall(i=1:mx)
                q(i,3) = q(i,3) + beta * exp(-((xlower+(i-0.5d0)*dx-location)/sigma)**2)
            end forall
            
!             ! Add perturbation of eigenvector to initial condition
!             ! Currently for 1st and 4th wave families
!             beta = 0.01d0
!             eig_vec(1) = 1.d0
!             one_minus_r = 1.d0 - r
!             do i=1,mx
!                 if (xlower+(i-0.5d0)*dx < -0.25d0) then
!                     mult_depth = q(i,1) * q(i,3)
! !                     lambda = - sqrt(g*total_depth) + 0.5d0 * mult_depth/ total_depth**(3/2) * one_minus_r
!                     lambda = sqrt(g*total_depth) + 0.5d0 * mult_depth/ total_depth**(3/2) * one_minus_r
!                     eig_vec(2) = lambda
!                     eig_vec(3) = (lambda**2 - g*q(i,1)) / (g*r*q(i,1))
!                     eig_vec(4) = eig_vec(3) * eig_vec(2)
!                     q(i,:) = q(i,:) + beta * eig_vec
!                 endif
!             enddo
    end select
             
    
end subroutine qinit