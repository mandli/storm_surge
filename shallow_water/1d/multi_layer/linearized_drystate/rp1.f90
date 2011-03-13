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
    double precision :: eigen_vectors(4,4),A(4,4),delta(4),alpha(4),beta(4)
    double precision, dimension(4) :: flux_l,flux_r
    double precision, dimension(2) :: h_hat_l,h_hat_r,h_tilde_l,h_tilde_r
    double precision, dimension(2) :: hu_tilde_l,hu_tilde_r
    double precision, dimension(2) :: h_l,u_l,hu_l,h_r,u_r,hu_r
    double precision :: b_l,b_r,gamma_l,gamma_r,h_ave(2)
    logical :: dry_state_l(2), dry_state_r(2)
    
    logical, parameter :: fwaves = .true.

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
        dry_state_l = .false.
        dry_state_r = .false.
        
        ! These state variables are the actual depth and momenta, note that
        ! the auxillary arrays contain the steady state depths that may be 
        ! spatial varying due to the bathymetry
        do j=1,2
            ! Steady states
            h_hat_l(j) = auxr(i-1,j+1)
            h_hat_r(j) = auxl(i,j+1)
            
            ! Perturbations
            h_tilde_l(j) = qr(i-1,(j-1)*2 + 1)
            hu_tilde_l(j) = qr(i-1,(j-1)*2 + 2)
            h_tilde_r(j) = ql(i,(j-1)*2 + 1)
            hu_tilde_r(j) = ql(i,(j-1)*2 + 2)
            
            ! Full variables
            h_l(j) = h_hat_l(j) + h_tilde_l(j)
            h_r(j) = h_hat_r(j) + h_tilde_r(j)
            hu_l(j) = hu_tilde_l(j)
            hu_r(j) = hu_tilde_r(j)
            
            ! Check for dry states in this layer
            if (h_l(j) < dry_tolerance) then
                dry_state_l(j) = .true.
            endif
            if (h_r(j) < dry_tolerance) then
                dry_state_r(j) = .true.
            endif
        enddo
        
        b_l = auxr(i-1,1)
        b_r = auxl(i,1)
        
        ! ====================================================================
        ! Calculate eigen-space based on initial steady state values in aux
        ! array
        !
        ! dry states for depth:
        !   h_l(1) = auxr(i-1,2)
        !   h_l(2) = auxr(i-1,3)
        !   h_r(1) = auxl(i,2)
        !   h_r(2) = auxl(i,3)
        ! These should be valid since we are assuming that the top layer does
        ! not go to zero
        gamma_l = h_hat_l(2) / h_hat_l(1)
        gamma_r = h_hat_r(2) / h_hat_r(1)
!         gamma_l = auxr(i-1,3) / auxr(i-1,2)
!         gamma_r = auxl(i,3) / auxl(i,2)
        
        alpha(1) = 0.5d0*(gamma_l-1+sqrt((gamma_l-1)**2+4.d0*r*gamma_l))
        alpha(2) = 0.5d0*(gamma_l-1-sqrt((gamma_l-1)**2+4.d0*r*gamma_l))
        alpha(3) = 0.5d0*(gamma_r-1-sqrt((gamma_r-1)**2+4.d0*r*gamma_r))
        alpha(4) = 0.5d0*(gamma_r-1+sqrt((gamma_r-1)**2+4.d0*r*gamma_r))

        s(i,1) = -sqrt(g*h_hat_l(1)*(1+alpha(1)))
        s(i,2) = -sqrt(g*h_hat_l(1)*(1+alpha(2)))
        s(i,3) = sqrt(g*h_hat_r(1)*(1+alpha(3)))
        s(i,4) = sqrt(g*h_hat_r(1)*(1+alpha(4)))
        
!         s(i,1) = -sqrt(g*auxr(i-1,2)*(1+alpha(1)))
!         s(i,2) = -sqrt(g*auxr(i-1,2)*(1+alpha(2)))
!         s(i,3) = sqrt(g*auxl(i,2)*(1+alpha(3)))
!         s(i,4) = sqrt(g*auxl(i,2)*(1+alpha(4)))
        
        eigen_vectors(1,:) = 1.d0
        eigen_vectors(2,:) = s(i,:)
        eigen_vectors(3,:) = alpha
        eigen_vectors(4,:) = s(i,:)*alpha
        
        if (dry_state_r(2).and.(.not.dry_state_l(2))) then
        endif
        
        ! ====================================================================
        ! Calculate flux vector to be projected onto e-space
        if (fwaves) then
            flux_r(1) = hu_r(1)
            flux_r(2) = hu_r(1)**2 / h_r(1) + 0.5d0 * g * h_r(1)**2
            flux_r(3) = hu_r(2)
            flux_r(4) = hu_r(2)**2 / h_r(2) + 0.5d0 * g * h_r(2)**2
            
            flux_l(1) = hu_l(1)
            flux_l(2) = hu_l(1)**2 / h_l(1) + 0.5d0 * g * h_l(1)**2
            flux_l(3) = hu_l(2)
            flux_l(4) = hu_l(2)**2 / h_l(2) + 0.5d0 * g * h_l(2)**2
            
            delta = flux_r - flux_l
            
            delta(2) = delta(2) + g * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
            delta(4) = delta(4) + g * h_ave(2) * (r * (h_r(1) - h_l(1)) + b_r - b_l)
            
!             ! [hu_1]
!             delta(1) = hu_tilde_r(1) - hu_tilde_l(1)
!             ! g ave(h_1) [h_1 + h_2]
!             delta(2) = 0.5d0 * g * (h_hat_r(1) + h_hat_l(1)) * (h_tilde_r(1) - h_tilde_l(1) + h_tilde_r(2) - h_tilde_l(2))
!             ! [hu_2]
!             delta(3) = hu_tilde_r(2) - hu_tilde_l(2)
!             ! g ave(h_2) ([h_2] + r [h_1])
!             
!             if (dry_state_r(2).and.(.not.dry_state_l(2))) then
!                 delta(4) = 0.5d0 * g * (h_tilde_r(2) + h_tilde_l(2)) * ((1.d0-r) * (b_r - b_l) + r * (h_tilde_r(1) - h_tilde_l(1)))
!             else
!                 delta(4) = 0.5d0 * g * (h_hat_r(2) + h_hat_l(2)) * (h_tilde_r(2) - h_tilde_l(2) + r * (h_tilde_r(1) - h_tilde_l(1)))
!             endif
            
            ! Add extra terms if near one sided dry state
!             if (.not.(dry_state_r(2).and.dry_state_l(2))) then
!                 if (dry_state_r(2).and.(.not.dry_state_l(2))) then
!                     eta_2 = h_hat_l(2) + b_l
!                     delta(2) = delta(2) - g * h_r(1) * (b_r - eta_2)
!                     delta(4) = delta(4) - r * g * h_l(2) * (b_r - eta_2)
!                 endif
!             endif

            !
            
            ! Wind forcing
!             wind_speed = 0.5d0 * (w_l + w_r)
!             tau = wind_drag(wind_speed) * rho_air * wind_speed
!             delta(4) = delta(4) - tau * wind_speed
        else
            delta = ql(i,:) - qr(i-1,:)
        endif


        ! ====================================================================
        ! Solve system, solution is stored in delta 
        print *,dry_state_l(2),dry_state_r(2)
        if (dry_state_l(2).and.dry_state_r(2)) then
            beta(1) = (delta(1)*s(i,4)-delta(2)) / (s(i,4)-s(i,1))
            beta(2) = 0.d0
            beta(3) = 0.d0
            beta(4) = (delta(2) - s(i,1)*beta(1)) / s(i,4)
        else if (dry_state_r(1).or.dry_state_l(1)) then
            stop "Top layer dry state not handled!  Exiting..."
        else
!             if (dry_state_r(2).and.(.not.dry_state_l(2))) then
!                 print *,"Here"
!             endif
            A = eigen_vectors
            call dgesv(4,1,A,4,ipiv,delta,4,info)
            if (.not.(info == 0)) then 
                print *, "Error solving R beta = delta, ",info
                print *, "Eigen-speeds:",(s(i,mw),mw=1,mwaves)
                print *, "Eigen-vectors:"
                do j=1,4
                    print "(4d16.8)",(eigen_vectors(j,m),m=1,meqn)
                enddo
                stop
            endif
            beta = delta
        endif
        
!         if (dry_state_r(2).and.(.not.dry_state_l(2))) then
!             print "('beta(',i3,')=',4d16.8)",i,beta
!             print *, "Eigen-speeds:",(s(i,mw),mw=1,mwaves)
!             print *, "Eigen-vectors:"
!             do j=1,4
!                 print "(4d16.8)",(eigen_vectors(j,m),m=1,meqn)
!             enddo
!         endif
        
        ! Calculate waves
        forall(mw=1:4)
            fwave(i,:,mw) = eigen_vectors(:,mw) * beta(mw)
        end forall

        ! Calculate amdq and apdq
        if (fwaves) then
            do mw=1,mwaves
                if (s(i,mw) > 0.d0) then
                    apdq(i,:) = apdq(i,:) + fwave(i,:,mw)
                else                                     
                    amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
                endif
            enddo
        else
            do mw=1,mwaves
                if (s(i,mw) > 0.d0) then
                    apdq(i,:) = apdq(i,:) + fwave(i,:,mw) * s(i,mw)
                else                                      
                    amdq(i,:) = amdq(i,:) + fwave(i,:,mw) * s(i,mw)
                endif
            enddo
        endif
    enddo

end subroutine rp1