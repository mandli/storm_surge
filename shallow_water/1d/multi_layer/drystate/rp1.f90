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
    double precision, dimension(2) :: h_l,u_l,hu_l,h_r,u_r,hu_r,h_ave(2)
    double precision :: b_l,b_r,gamma_l,gamma_r,momentum_transfer
    logical :: dry_state_l(2), dry_state_r(2)
    
    double precision :: total_depth_l,total_depth_r,mult_depth_l,mult_depth_r 
    double precision :: kappa_l,kappa_r

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
        
        ! These state variables are the actual depth and momenta
        do j=1,2
            h_l(j) = qr(i-1,2*j-1)
            h_r(j) = ql(i  ,2*j-1)
            hu_l(j) = qr(i-1,2*j)
            hu_r(j) = ql(i  ,2*j)
            
            ! Check for dry states in this layer
            if (h_l(j) < dry_tolerance) then
                dry_state_l(j) = .true.
                u_l(j) = 0.d0
            else
                u_l(j) = hu_l(j) / h_l(j)
            endif
            if (h_r(j) < dry_tolerance) then
                dry_state_r(j) = .true.
                u_r(j) = 0.d0
            else
                u_r(j) = hu_r(j) / h_r(j)
            endif
        enddo
        
        b_l = auxr(i-1,1)
        b_r = auxl(i,1)
        
        ! ====================================================================
        ! Calculate eigen-space values
        gamma_l = h_l(2) / h_l(1)
        gamma_r = h_r(2) / h_r(1)
        
        total_depth_l = sum(h_l)
        total_depth_r = sum(h_r)
        mult_depth_l = product(h_l)
        mult_depth_r = product(h_r)
        
        kappa_l = (u_l(1)-u_l(2))**2 / (g*one_minus_r*total_depth_l)
        kappa_r = (u_r(1)-u_r(2))**2 / (g*one_minus_r*total_depth_r)
        if ((kappa_l > 1.d0).or.(kappa_r > 1.d0)) then
            print "(a,2d16.8)","Hyperbolicity may have failed, kappa = ",kappa_l,kappa_r
        endif
        
        alpha(1) = 0.5d0*(gamma_l-1.d0+sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
        alpha(2) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
        alpha(3) = 0.5d0*(gamma_r-1.d0-sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
        alpha(4) = 0.5d0*(gamma_r-1.d0+sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))

        s(i,1) = -sqrt(g*h_l(1)*(1+alpha(1)))
        s(i,2) = -sqrt(g*h_l(1)*(1+alpha(2)))
        s(i,3) = sqrt(g*h_r(1)*(1+alpha(3)))
        s(i,4) = sqrt(g*h_r(1)*(1+alpha(4)))

!         s(i,1) = (h_l(1)*u_l(1)+h_l(2)*u_l(2)) / total_depth_l - sqrt(g*total_depth_l)
!         s(i,2) = (h_l(2)*u_l(1)+h_l(1)*u_l(2)) / total_depth_l - sqrt(g*one_minus_r*mult_depth_l/total_depth_l * (1-(u_l(1)-u_l(2))**2/(g*one_minus_r*total_depth_l)))
!         s(i,3) = (h_r(2)*u_r(1)+h_r(1)*u_r(2)) / total_depth_r + sqrt(g*one_minus_r*mult_depth_r/total_depth_r * (1-(u_r(1)-u_r(2))**2/(g*one_minus_r*total_depth_r)))
!         s(i,4) = (h_r(1)*u_r(1)+h_r(2)*u_r(2)) / total_depth_r + sqrt(g*total_depth_r)
        
!         alpha(1:2) = g * h_l(2) / (s(i,1:2)**2 - g*h_l(2))
!         alpha(3:4) = g * h_r(2) / (s(i,3:4)**2 - g*h_r(2))
        
        eigen_vectors(1,:) = 1.d0
        eigen_vectors(2,:) = s(i,:)
        eigen_vectors(3,:) = alpha
        eigen_vectors(4,:) = s(i,:)*alpha(:)

        
        
        ! ====================================================================
        ! Calculate flux vector to be projected onto e-space
        if (fwaves) then
            do j=1,2
                flux_r(2*j-1) = hu_r(j)
                flux_l(2*j-1) = hu_l(j)
                flux_r(2*j) = h_r(j) * u_r(j)**2 + 0.5d0 * g * h_r(j)**2
                flux_l(2*j) = h_l(j) * u_l(j)**2 + 0.5d0 * g * h_l(j)**2
            enddo

            delta = flux_r - flux_l
            
            h_ave = 0.5d0 * (h_r + h_l)
            delta(2) = delta(2) + g * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
            delta(4) = delta(4) + g * h_ave(2) * (r * (h_r(1) - h_l(1)) + b_r - b_l)

            ! Wind forcing
!             wind_speed = 0.5d0 * (w_l + w_r)
!             tau = wind_drag(wind_speed) * rho_air * wind_speed
!             delta(4) = delta(4) - tau * wind_speed
        else
            delta = ql(i,:) - qr(i-1,:)
        endif


        ! ====================================================================
        ! Solve system, solution is stored in delta
        ! Should we check here that delta(4) is zero?
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