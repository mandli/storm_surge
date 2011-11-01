! ============================================================================
!  Reimann Solver
! ============================================================================
subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
!
!   Riemann solver for linearized multilayer shallow water equations
!

    use parameters_module

    implicit none
    
    ! Input arguments
    integer, intent(in) :: maxmx,meqn,mwaves,mbc,mx
    
    double precision, intent(in), dimension(1-mbc:maxmx+mbc, meqn) :: ql,qr
    double precision, intent(in), dimension(1-mbc:maxmx+mbc, *) :: auxl,auxr
    
    ! Output arguments
    double precision, intent(out) :: s(1-mbc:maxmx+mbc, mwaves)
    double precision, intent(out) :: fwave(1-mbc:maxmx+mbc, meqn, mwaves)
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn) :: amdq,apdq
    
    ! Locals
    integer :: i,j,m,mw,ipiv(4),info
    double precision :: eig_vec(4,4),A(4,4),delta(4),alpha(4),beta(4)
    double precision, dimension(4) :: flux_l,flux_r
    double precision, dimension(2) :: h_l,u_l,hu_l,h_r,u_r,hu_r,h_ave,u_ave
    double precision :: b_l,b_r,gamma_l,gamma_r,tau,w_l,w_r
    double precision :: wind_speed,lambda(4),eta_l(2),eta_r(2),h_hat_l(2),h_hat_r(2)
    logical :: dry_state_l(2), dry_state_r(2), inundation
    
    integer :: layer_index,rare(1-mbc:mx+mbc)
    double precision :: momentum_transfer(2),flux_transfer_r,flux_transfer_l
    double precision :: inundation_height(2)

    ! Eigenvalue functions
    integer :: single_layer_eigen,linear_eigen,velocity_eigen,lapack_eigen

    ! Common block
    double precision :: dt,dx,t
    common /comxt/ dt,dx,t

    ! Initialize return variables
    amdq = 0.d0
    apdq = 0.d0

    !          |          |          
    !          |          |          
    !          |          |          
    !----------|----------|----------
    !    i-1         i         i+1
    !  qr(i-1)  ql(i)
    !     Riemann problem at i-1/2 is between qr(i-1) (left state) and ql(i) 
    !     (right state)
    do i=2-mbc,mx+mbc
        inundation = .false.
        dry_state_l = .false.
        dry_state_r = .false.
        
        ! These state variables are the actual depth and momenta
        do j=1,2
            layer_index = 2*(j-1)
            h_l(j) = qr(i-1,layer_index+1) / rho(j)
            h_r(j) = ql(i,layer_index+1) / rho(j)
            hu_l(j) = qr(i-1,layer_index+2) / rho(j)
            hu_r(j) = ql(i,layer_index+2) / rho(j)
            
            h_hat_l(j) = auxr(i-1,j+2)
            h_hat_r(j) = auxl(i,j+2)
            
            ! Check for dry states in this layer
            if (h_l(j) < dry_tolerance) then
                dry_state_l(j) = .true.
                h_l(j) = 0.d0
                u_l(j) = 0.d0
            else
                u_l(j) = qr(i-1,layer_index+2) / qr(i-1,layer_index+1)
            endif
            if (h_r(j) < dry_tolerance) then
                dry_state_r(j) = .true.
                h_r(j) = 0.d0
                u_r(j) = 0.d0
            else
                u_r(j) = ql(i,layer_index+2) / ql(i,layer_index+1)
            endif
        enddo
        
        h_ave = 0.5d0 * (h_l(:) + h_r(:))
        
        b_l = auxr(i-1,1)
        b_r = auxl(i,1)
        w_l = auxr(i-1,2)
        w_r = auxl(i,2)      
            
        ! ====================================================================
        ! Solve Single layer problem seperately
        if (dry_state_r(2).and.dry_state_l(2)) then
            rare(i) = single_layer_eigen(h_l,h_r,u_l,u_r,b_l,b_r,lambda,eig_vec)
            s(i,:) = lambda
            
            delta(1) = rho(1) * (hu_r(1) - hu_l(1))
            flux_r(2) = rho(1) * (h_r(1) * u_r(1)**2 + 0.5d0 * g * h_r(1)**2)
            flux_l(2) = rho(1) * (h_l(1) * u_l(1)**2 + 0.5d0 * g * h_l(1)**2)
            delta(2) = flux_r(2) - flux_l(2) + g * rho(1) * h_ave(1) * (b_r - b_l)
            
            beta(1) = (delta(1)*s(i,4)-delta(2)) / (s(i,4)-s(i,1))
            beta(2) = 0.d0
            beta(3) = 0.d0
            beta(4) = (delta(2) - s(i,1)*beta(1)) / s(i,4)
            
            ! Calculate waves
            forall(mw=1:4)
                fwave(i,:,mw) = eig_vec(:,mw) * beta(mw)
            end forall
            cycle
        endif
        
        ! ====================================================================
        ! Calculate eigen-space values
        ! ====================================================================
        
        ! ====================================================================
        ! Inundation cases
        if (dry_state_r(2).and.(.not.dry_state_l(2)).and.(h_l(2) + b_l > b_r)) then
            inundation = .true.
            print *,"Right inundation problem"
            inundation = .true.
            if (inundation_method == 0) then
                stop "Inundation not allowed."
            else if (inundation_method == 1) then
                ! Linear eigensystem
                inundation_height = [h_r(1),0.d0]
                rare(i) = linear_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
                
                ! Corrective wave
                s(i,3) = u_l(2) + 2.d0 * sqrt(g*(1.d0-r)*h_l(2))
                s(i,4) = u_r(1) + sqrt(g*h_r(1))
                alpha(3) = r * g * h_l(2) / ((s(i,3) - u_l(2))**2 - g * h_l(2))
                alpha(4) = 0.d0

                eig_vec(1,3:4) = 1.d0
                eig_vec(2,3:4) = s(i,3:4)
                eig_vec(3,3:4) = alpha(3:4)
                eig_vec(4,3:4) = s(i,3:4)*alpha(3:4)
            else if (inundation_method == 2) then
                inundation_height = [h_r(1),dry_tolerance]
                h_r(2) = dry_tolerance
                rare(i) = linear_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
                
                s(i,4) = u_r(1) + sqrt(g*h_r(1))
                eig_vec(:,4) = [1.d0,s(i,4),0.d0,0.d0]
            else if (inundation_method == 3) then
                inundation_height = [h_r(1),dry_tolerance]
                rare(i) = velocity_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
                
                ! Correction for fast wave
                s(i,4) = u_r(1) + sqrt(g*h_r(1))
                eig_vec(:,4) = [1.d0,s(i,4),0.d0,0.d0]
            else if (inundation_method == 4) then
                inundation_height = [h_r(1),dry_tolerance]
                rare(i) = lapack_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
                
                ! Correction for the fast waves
                s(i,2) = u_r(1) + sqrt(g*h_r(1))
                eig_vec(:,2) = [1.d0,s(i,2),0.d0,0.d0]
            else if (inundation_method == 5) then
                inundation_height = [h_r(1),0.d0]
                rare(i) = lapack_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda                
            endif   
        else if (dry_state_l(2).and.(.not.dry_state_r(2)).and.(h_r(2) + b_r > b_l)) then
            inundation = .true.
            print *,"Left inundation problem"
            inundation = .true.
            ! Inundation problem eigen
            if (inundation_method == 0) then
                stop "Inundation not allowed."
            else if (inundation_method == 1) then
                ! Linear eigensystem
                inundation_height = [h_l(1),dry_tolerance]
                rare(i) = linear_eigen(inundation_height,h_r,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda

                ! Corrections to internal wave
                s(i,2) = u_r(2) - 2.d0 * sqrt(g*(1.d0-r)*h_r(2))
                alpha(2) = r * g * h_r(2) / ((s(i,2) - u_r(2))**2 - g * h_r(2))
                eig_vec(:,2) = [1.d0,s(i,2),alpha(2),alpha(2)*s(i,2)]
                
                ! Correction for the fast waves
                s(i,1) = u_l(1) - sqrt(g*h_l(1))
                eig_vec(:,1) = [1.d0,s(i,1),0.d0,0.d0]
            else if (inundation_method == 2) then
                ! Use linearized eigensystem
                inundation_height = [h_l(1),dry_tolerance]
                rare(i) = linear_eigen(inundation_height,h_r,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
                            
                ! Correction for the fast waves
                s(i,1) = u_l(1) - sqrt(g*h_l(1))
                eig_vec(:,1) = [1.d0,s(i,1),0.d0,0.d0]
            else if (inundation_method == 3) then
                ! Use velocity difference expansion eigensystems
                inundation_height = [h_l(1),dry_tolerance]
                rare(i) = velocity_eigen(inundation_height,h_r,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda                
        
                ! Correction for the fast waves
                s(i,1) = u_r(1) - sqrt(g*h_r(1))
                eig_vec(:,1) = [1.d0,s(i,1),0.d0,0.d0]
            else if (inundation_method == 4) then
                ! LAPACK solver with corrective wave and small wet layer
                inundation_height = [h_r(1),dry_tolerance]
                rare(i) = lapack_eigen(h_l,inundation_height,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
                            
                ! Correction for the fast waves
                s(i,1) = u_l(1) - sqrt(g*h_l(1))
                eig_vec(:,1) = [1.d0,s(i,1),0.d0,0.d0]
            else if (inundation_method == 5) then
                ! Use the LAPACK solver with no correction
                inundation_height = [h_l(1),0.d0]
                rare(i) = lapack_eigen(inundation_height,h_r,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
            endif        
          
        ! ====================================================================  
        !  Wall or wet case       
        else
            ! Wall dry state or completely wet case
            if (eigen_method == 1) then
                rare(i) = linear_eigen(h_hat_l,h_hat_r,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
            else if (eigen_method == 2) then
                rare(i) = linear_eigen(h_l,h_r,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
            else if (eigen_method == 3) then 
                rare(i) = velocity_eigen(h_l,h_r,u_l,u_r,b_l,b_r,lambda,eig_vec)
                s(i,:) = lambda
            else if (eigen_method == 4) then
                if (dry_state_r(2).and.(.not.dry_state_l(2)).or. &
                        dry_state_l(2).and.(.not.dry_state_r(2))) then
                    rare(i) = linear_eigen(h_l,h_r,u_l,u_r,b_l,b_r,lambda,eig_vec)
                    s(i,:) = lambda
                else
                    rare(i) = lapack_eigen(h_l,h_r,u_l,u_r,b_l,b_r,lambda,eig_vec)
                    s(i,:) = lambda
                endif
            else
                stop "Invalid eigensystem method requested, method = (1,4)."
            endif
        endif
        
        !  end of eigenspace calculation
        ! ====================================================================
        
        ! ====================================================================
        ! Calculate flux vector to be projected onto e-space
        ! Calculate jumps in fluxes
        if (dry_state_r(2).and.(.not.dry_state_l(2)).and.(.not.inundation)) then
            ! Wall boundary conditions
            h_r(2) = h_l(2)
            hu_r(2) = -hu_l(2)
            u_r(2) = -u_l(2)
            
            ! Top layer eta(2) = b_r - h_l(2) - b_l
            momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r - h_l(2) - b_l)
            
            momentum_transfer(2) = 0.d0
            flux_transfer_r = 0.d0
            flux_transfer_l = 0.d0
        ! Left state dry, right wet
        else if (dry_state_l(2).and.(.not.dry_state_r(2)).and.(.not.inundation)) then
            ! Wall boundary conditions
            h_l(2) = h_r(2)
            hu_l(2) = -hu_r(2)
            u_l(2) = -u_r(2)
            
            ! Top layer eta(2) = h_r(2) + b_r - b_l
            momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r + h_r(2) - b_l)
            
            momentum_transfer(2) = 0.d0
            flux_transfer_r = 0.d0
            flux_transfer_l = 0.d0
        ! Fully wet bottom layer or inundation
        else
            momentum_transfer(1) =   g * rho(1) * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
            momentum_transfer(2) = - g * rho(1) * h_ave(1) * (h_r(2) - h_l(2)) + g * rho(2) * h_ave(2) * (b_r - b_l)
            ! Bottom layer momentum transfer flux
            flux_transfer_r = g * rho(1) * product(h_r)
            flux_transfer_l = g * rho(1) * product(h_l)
        endif
        
        ! Flux jumps
        do j=1,2
            layer_index = 2*(j-1)
            flux_r(layer_index+1) = rho(j) * hu_r(j)
            flux_r(layer_index+2) = rho(j) * (h_r(j) * u_r(j)**2 + 0.5d0 * g * h_r(j)**2)
            
            flux_l(layer_index+1) = rho(j) * hu_l(j)
            flux_l(layer_index+2) = rho(j) * (h_l(j) * u_l(j)**2 + 0.5d0 * g * h_l(j)**2)
        enddo
        
        delta = flux_r - flux_l

        ! Bottom layer additional flux
        delta(4) = delta(4) + flux_transfer_r - flux_transfer_l

        ! Momentum source term from layers
        delta(2) = delta(2) + momentum_transfer(1)
        delta(4) = delta(4) + momentum_transfer(2)
        
        ! Wind forcing
        wind_speed = 0.5d0 * (w_l + w_r)
        tau = wind_drag(wind_speed) * rho_air * wind_speed
        delta(2) = delta(2) - tau * wind_speed

        ! ====================================================================
        ! Solve system, solution is stored in delta
        A = eig_vec
        call dgesv(4,1,A,4,ipiv,delta,4,info)
        if (.not.(info == 0)) then 
            print *, "Location (i) = (",i,")"
            print *, "Dry states, L=",dry_state_l(2)," R=",dry_state_r(2)
            print *, "h_l(2) = ",h_l(2)," h_r(2) = ",h_r(2)
            print *, "Error solving R beta = delta, ",info
            print *, "Eigen-speeds:",(s(i,mw),mw=1,mwaves)
            print *, "Eigen-vectors:"
            do j=1,4
                print "(4d16.8)",(eig_vec(j,m),m=1,meqn)
            enddo
            stop
        endif
        beta = delta
        
        ! Calculate waves
        forall(mw=1:4)
            fwave(i,:,mw) = eig_vec(:,mw) * beta(mw)
        end forall
    enddo
    
    ! Calculate amdq and apdq
    ! No entropy fix requested
    if (.not.entropy_fix) then
        do i=2-mbc,mx+mbc
            do mw=1,mwaves
                if (s(i,mw) > 0.d0) then
                    apdq(i,:) = apdq(i,:) + fwave(i,:,mw)
                else                                     
                    amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
                endif
            enddo
        enddo
    else
        ! if rare(i) > 0.d0 then wave correction is in 2nd family
        ! if rare(i) < 0.d0 then wave correction is in 3rd family
        do i=2-mbc,mx+mbc
            ! Move from left to right acros wave fan checking position of each
            ! wave, here we only need to figure out the left going
            ! fluctuations amdq since we can calculate apdq from that quantity
            ! as     apdq = sum fwave - amdq     holds
            
            ! (1) Fully supersonic case right going 
            if (s(i,1) > 0.d0) then
                continue
            ! (2) All but the 1st wave is going right
            else if (s(i,2) > 0.d0) then
                amdq(i,:) = amdq(i,:) + fwave(i,:,1)
            ! (3) Rarefaction in 2nd wave family
            else if (rare(i) > 0.d0) then
                ! Note that this will be zero unless this is an inundation
                ! problem, puzzling
                beta = (lambda_r - s(i,2)) / (lambda_r - lambda_l)
                amdq(i,:) = amdq(i,:) + fwave(i,:,1)
                amdq(i,:) = amdq(i,:) + beta*fwave(i,:,2)
            ! (4) 1st and 2nd waves going right
            else if (s(i,2) < 0.d0) then
                do mw=1,2
                    amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
                enddo
            ! (5) Rarefaction in 3rd family
            else if (rare(i) < 0.d0) then
                do mw=1,2
                    amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
                enddo
                lambda_l = rare(i)
                lambda_r = s(i,3)
                beta = (lambda_r - s(i,3)) / (lambda_r - lambda_l)
                amdq(i,:) = amdq(i,:) + beta * fwave(i,:,3)
            ! (6) 1st, 2nd, and 3rd going right
            else if (s(i,3) < 0.d0) then
                do mw=1,3
                    amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
                enddo
            ! (7) Fully supersonic case left going
            else if (s(i,4) < 0.d0) then
                do mw=1,mwaves
                    amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
                enddo
            endif
            
            ! Compute each element going right
            do m=1,meqn
                apdq(i,m) = sum(fwave(i,m,:)) - amdq(i,m)
            enddo
        enddo
    endif
end subroutine rp1


! ============================================================================
!  Eigenvalue routines
! ============================================================================
double precision function linear_eigen(h_l,h_r,u_l,u_r,b_l,b_r,s,eig_vec) result(rare)

    use parameters_module, only: r,g,dry_tolerance,entropy_fix

    implicit none
    
    ! I/O
    double precision, intent(in) :: h_l(2),h_r(2),u_l(2),u_r(2),b_l,b_r
    double precision, intent(inout) :: s(4),eig_vec(4,4)
    
    ! Locals
    double precision :: alpha(4),gamma_l,gamma_r
        
    rare = 0.d0
        
    gamma_l = h_l(2) / h_l(1)
    gamma_r = h_r(2) / h_r(1)

    alpha(1) = 0.5d0*(gamma_l-1.d0+sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
    alpha(2) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
    alpha(3) = 0.5d0*(gamma_r-1.d0-sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
    alpha(4) = 0.5d0*(gamma_r-1.d0+sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))

    s(1) = u_l(1) - sqrt(g*h_l(1)*(1+alpha(1)))
    s(2) = u_l(2) - sqrt(g*h_l(1)*(1+alpha(2)))
    s(3) = u_r(2) + sqrt(g*h_r(1)*(1+alpha(3)))
    s(4) = u_r(2) + sqrt(g*h_r(1)*(1+alpha(4)))

    eig_vec(1,:) = 1.d0
    eig_vec(2,:) = s(:)
    eig_vec(3,:) = alpha
    eig_vec(4,:) = s(:)*alpha(:)

    if (entropy_fix) then
        ! Check for entropy problems, reverse state evaluation
        alpha(2) = 0.5d0*(gamma_r-1.d0-sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
        alpha(3) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
    
        if (u_r(2) - sqrt(g*h_r(1)*(1+alpha(2))) > 0.d0) then
            ! Entropy correction required for wave 2
            rare = u_r(2) - sqrt(g*h_r(1)*(1+alpha(2)))
        else if (u_l(2) + sqrt(g*h_l(1)*(1+alpha(3))) < 0.d0) then
            ! Entropy correction required for wave 3
            rare = u_l(2) + sqrt(g*h_l(1)*(1+alpha(3)))
        endif
    endif

end function linear_eigen

double precision function velocity_eigen(h_l,h_r,u_l,u_r,b_l,b_r,s,eig_vec) result(rare)

    use parameters_module, only: r,g,one_minus_r,entropy_fix

    implicit none
    
    ! I/O
    double precision, intent(in) :: h_l(2),h_r(2),u_l(2),u_r(2),b_l,b_r
    double precision, intent(inout) :: s(4),eig_vec(4,4)
    
    ! Locals
    double precision :: total_depth_l,total_depth_r,mult_depth_l,mult_depth_r
    double precision :: alpha(4)
    
    rare = 0
    
    total_depth_l = sum(h_l)
    total_depth_r = sum(h_r)
    mult_depth_l = product(h_l)
    mult_depth_r = product(h_r)

    s(1) = (h_l(1)*u_l(1)+h_l(2)*u_l(2)) / total_depth_l &
          - sqrt(g*total_depth_l)
    s(2) = (h_l(2)*u_l(1)+h_l(1)*u_l(2)) / total_depth_l &
          - sqrt(g*one_minus_r*mult_depth_l/total_depth_l &
          * (1-(u_l(1)-u_l(2))**2/(g*one_minus_r*total_depth_l)))
    s(3) = (h_r(2)*u_r(1)+h_r(1)*u_r(2)) / total_depth_r &
          + sqrt(g*one_minus_r*mult_depth_r/total_depth_r &
          * (1-(u_r(1)-u_r(2))**2/(g*one_minus_r*total_depth_r)))
    s(4) = (h_r(1)*u_r(1)+h_r(2)*u_r(2)) / total_depth_r &
            + sqrt(g*total_depth_r)

    alpha(1:2) = ((s(1:2) - u_l(1))**2 - g * h_l(1)) / (g*h_l(1))
    alpha(3:4) = ((s(3:4) - u_r(1))**2 - g * h_r(1)) / (g*h_r(1))

    eig_vec(1,:) = 1.d0
    eig_vec(2,:) = s(:)
    eig_vec(3,:) = alpha
    eig_vec(4,:) = s(:)*alpha(:)
    

    if (entropy_fix) then
        ! Check for entropy problems, reverse state evaluation
        alpha(2) = (h_r(2)*u_r(1)+h_r(1)*u_r(2)) / total_depth_r &
              - sqrt(g*one_minus_r*mult_depth_r/total_depth_r &
              * (1-(u_r(1)-u_r(2))**2/(g*one_minus_r*total_depth_r)))
        alpha(3) = (h_l(2)*u_l(1)+h_l(1)*u_l(2)) / total_depth_l &
              + sqrt(g*one_minus_r*mult_depth_l/total_depth_l &
              * (1-(u_l(1)-u_l(2))**2/(g*one_minus_r*total_depth_l)))
        if (alpha(2) > 0.d0) then
            rare = alpha(2)
        else if (alpha(3) < 0.d0) then
            rare = alpha(3)
        endif
    endif

end function velocity_eigen


double precision function lapack_eigen(h_l,h_r,u_l,u_r,b_l,b_r,s,eig_vec) result(rare)

    use parameters_module, only: r,g,entropy_fix

    implicit none
    
    ! I/O
    double precision, intent(in) :: h_l(2),h_r(2),u_l(2),u_r(2),b_l,b_r
    double precision, intent(inout) :: s(4),eig_vec(4,4)
    
    ! Local
    integer, parameter :: lwork = 4*4
    integer :: j,info
    double precision :: A(4,4),h_ave(2),u_ave(2)
    double precision :: real_evalues(4),imag_evalues(4)
    double precision :: empty,work(1,lwork)
    
    ! Solve eigenvalue problem
    h_ave(:) = 0.5d0 * (h_l(:) + h_r(:))
    u_ave(:) = 0.5d0 * (u_l(:) + u_r(:))
    A(1,:) = [0.d0,1.d0,0.d0,0.d0]
    A(2,:) = [-u_ave(1)**2 + g*h_ave(1),2.d0*u_ave(1),g*h_ave(1),0.d0]
    A(3,:) = [0.d0,0.d0,0.d0,1.d0]
    A(4,:) = [g*r*h_ave(2),0.d0,-u_ave(2)**2 + g*h_ave(2),2.d0*u_ave(2)]
    call dgeev('N','V',4,A,4,real_evalues,imag_evalues,empty,1,eig_vec,4,work,lwork,info)
    if (info < 0) then
        info = -info
        print "(a,i1,a)","The ",info,"th argument had an illegal value."
        stop
    else if (info > 0) then
        print "(a)","The QR algorithm failed to compute all the"
        print "(a)","eigenvalues, and no eigenvectors have been"
        print "(a,i1,a)","computed; elements",info,"+1:4 of WR and WI"
        print "(a)","contain eigenvalues which have converged."
        stop
    endif
    do j=1,4
!         if (eig_vec(1,j) /= 0.d0) then
!             eig_vec(:,j) = eig_vec(:,j) / eig_vec(1,j)
!         endif
        if (imag_evalues(j) > 0.d0) then
            print "(a,i1,a,d16.8)","Imaginary eigenvalue(",j,") > 0.0",imag_evalues(j)
            stop
        endif
        s(j) = real_evalues(j)
    enddo
    
    rare = 0.d0
    if (entropy_fix) then
        stop "ERROR:  Entropy fix has not been implemented with this solver."
    endif
    
end function lapack_eigen

double precision function single_layer_eigen(h_l,h_r,u_l,u_r,b_l,b_r,s,eig_vec) result(rare)

    use parameters_module, only: g

    implicit none
    
    ! I/O
    double precision, intent(in) :: h_l(2),h_r(2),u_l(2),u_r(2),b_l,b_r
    double precision, intent(inout) :: s(4),eig_vec(4,4)
    
    s(1) = u_l(1) - sqrt(g*h_l(1))
    s(2) = 0.d0
    s(3) = 0.d0
    s(4) = u_r(1) + sqrt(g*h_r(1))
    
    eig_vec(1,:) = 1.d0
    eig_vec(2,:) = s(:)
    eig_vec(3,:) = 0.d0
    eig_vec(4,:) = 0.d0
    
    rare = 0.d0

end function single_layer_eigen
! ============================================================================

