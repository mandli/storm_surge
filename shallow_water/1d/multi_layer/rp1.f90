subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
!
!   Riemann solver for linearized multilayer shallow water equations
!

    use parameters_module

    implicit none
    
    ! Input arguments
    integer, intent(in) :: maxmx,meqn,mwaves,mbc,mx
    
    double precision, intent(in), dimension(1-mbc:maxmx+mbc, meqn) :: ql,qr
    double precision, intent(inout), dimension(1-mbc:maxmx+mbc, *) :: auxl,auxr
    
    ! Output arguments
    double precision, intent(out) :: s(1-mbc:maxmx+mbc, mwaves)
    double precision, intent(out) :: fwave(1-mbc:maxmx+mbc, meqn, mwaves)
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn) :: amdq,apdq
    
    ! Locals
    integer :: i,j,m,mw,ipiv(4),info
    double precision :: eigen_vectors(4,4),A(4,4),delta(4),alpha(4),beta(4)
    double precision, dimension(4) :: flux_l,flux_r
    double precision, dimension(2) :: h_l,u_l,hu_l,h_r,u_r,hu_r,h_ave,u_ave
    double precision :: b_l,b_r,gamma_l,gamma_r,momentum_transfer,tau,w_l,w_r
    double precision :: kappa_l,kappa_r,total_depth_l,total_depth_r,wind_speed
    double precision :: mult_depth_l,mult_depth_r,eta_l(2),eta_r(2)
    logical :: dry_state_l(2), dry_state_r(2)
        
    integer, parameter :: lwork = 4*4
    double precision :: real_evalues(4),imag_evalues(4)
    double precision :: empty,work(1,lwork)
    
    double precision, parameter :: RICH_TOLERANCE = 0.95d0

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
        dry_state_l = .false.
        dry_state_r = .false.
        
        ! These state variables are the actual depth and momenta
        do j=1,2
            h_l(j) = qr(i-1,2*j-1) / rho(j)
            h_r(j) = ql(i,2*j-1) / rho(j)
            hu_l(j) = qr(i-1,2*j) / rho(j)
            hu_r(j) = ql(i,2*j) / rho(j)
            
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
        w_l = auxr(i-1,2)
        w_r = auxl(i,2)
        
        kappa_l = (u_l(1)-u_l(2))**2 / (g*one_minus_r*sum(h_l))
        kappa_r = (u_r(1)-u_r(2))**2 / (g*one_minus_r*sum(h_r))
        if ((kappa_l > RICH_TOLERANCE).and.(.not.dry_state_l(2))) then
            print "(a,i4,a,d16.8)","Hyperbolicity may have failed, kappa(",i,") = ",kappa_l
        else if ((kappa_r > RICH_TOLERANCE).and.(.not.dry_state_r(2))) then
            print "(a,i4,a,d16.8)","Hyperbolicity may have failed, kappa(",i,") = ",kappa_r
        endif
        
        auxl(i,5) = kappa_l
        auxr(i-1,5) = kappa_r
        
        ! ====================================================================
        ! Calculate eigen-space values
        if ((eigen_method > 0).and.(eigen_method <= 3)) then
            if (eigen_method == 1) then
                gamma_l = auxr(i-1,4) / auxr(i-1,3)
                gamma_r = auxl(i,4) / auxl(i,3)
            
                alpha(1) = 0.5d0*(gamma_l-1.d0+sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
                alpha(2) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
                alpha(3) = 0.5d0*(gamma_r-1.d0-sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
                alpha(4) = 0.5d0*(gamma_r-1.d0+sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))

                s(i,1) = -sqrt(g*auxr(i-1,3)*(1+alpha(1)))
                s(i,2) = -sqrt(g*auxr(i-1,3)*(1+alpha(2)))
                s(i,3) = sqrt(g*auxl(i,3)*(1+alpha(3)))
                s(i,4) = sqrt(g*auxl(i,3)*(1+alpha(4)))
            else if (eigen_method == 2) then
                gamma_l = h_l(2) / h_l(1)
                gamma_r = h_r(2) / h_r(1)
            
                alpha(1) = 0.5d0*(gamma_l-1.d0+sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
                alpha(2) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
                alpha(3) = 0.5d0*(gamma_r-1.d0-sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
                alpha(4) = 0.5d0*(gamma_r-1.d0+sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))

                s(i,1) = -sqrt(g*h_l(1)*(1+alpha(1)))
                s(i,2) = -sqrt(g*h_l(1)*(1+alpha(2)))
                s(i,3) = sqrt(g*h_r(1)*(1+alpha(3)))
                s(i,4) = sqrt(g*h_r(1)*(1+alpha(4)))
            else if (eigen_method == 3) then 
                if (dry_state_l(2).and.dry_state_r(2)) then
                    s(i,1) = u_l(1) - sqrt(g*h_l(1))
                    s(i,2) = 0.d0
                    s(i,3) = 0.d0
                    s(i,4) = u_r(1) + sqrt(g*h_r(1))
                    alpha = 0.d0
                else
                    total_depth_l = sum(h_l)
                    total_depth_r = sum(h_r)
                    mult_depth_l = product(h_l)
                    mult_depth_r = product(h_r)

                    s(i,1) = (h_l(1)*u_l(1)+h_l(2)*u_l(2)) / total_depth_l &
                                - sqrt(g*total_depth_l)
                    s(i,2) = (h_l(2)*u_l(1)+h_l(1)*u_l(2)) / total_depth_l &
                                - sqrt(g*one_minus_r*mult_depth_l/total_depth_l &
                                * (1-(u_l(1)-u_l(2))**2/(g*one_minus_r*total_depth_l)))
                    s(i,3) = (h_r(2)*u_r(1)+h_r(1)*u_r(2)) / total_depth_r &
                                + sqrt(g*one_minus_r*mult_depth_r/total_depth_r &
                                * (1-(u_r(1)-u_r(2))**2/(g*one_minus_r*total_depth_r)))
                    s(i,4) = (h_r(1)*u_r(1)+h_r(2)*u_r(2)) / total_depth_r &
                                + sqrt(g*total_depth_r)
                
                    alpha(1:2) = ((s(i,1:2) - u_l(1))**2 - g * h_l(1)) / (g*h_l(1))
                    alpha(3:4) = ((s(i,3:4) - u_r(1))**2 - g * h_r(1)) / (g*h_r(1))
                endif
            endif
        
            eigen_vectors(1,:) = 1.d0
            eigen_vectors(2,:) = s(i,:)
            eigen_vectors(3,:) = alpha
            eigen_vectors(4,:) = s(i,:)*alpha(:)
        else if (eigen_method == 4) then
            if (dry_state_l(2).and.dry_state_r(2)) then
                s(i,1) = u_l(1) - sqrt(g*h_l(1))
                s(i,2) = 0.d0
                s(i,3) = 0.d0
                s(i,4) = u_r(1) + sqrt(g*h_r(1))
                eigen_vectors(1,:) = 1.d0
                eigen_vectors(2,:) = s(i,:)
                eigen_vectors(3,:) = 0.d0
                eigen_vectors(4,:) = 0.d0
            else if (dry_state_r(2).and.(.not.dry_state_l(2))) then
                gamma_l = auxr(i-1,4) / auxr(i-1,3)
                gamma_r = auxl(i,4) / auxl(i,3)
            
                alpha(1) = 0.5d0*(gamma_l-1.d0+sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
                alpha(2) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
                
                s(i,1) = -sqrt(g*auxr(i-1,3)*(1+alpha(1)))
                s(i,2) = -sqrt(g*auxr(i-1,3)*(1+alpha(2)))
                s(i,3) = 0.d0
                s(i,4) = u_r(1) + sqrt(g*h_r(1))
                eigen_vectors(1,1:2) = 1.d0
                eigen_vectors(2,1:2) = s(i,1:2)
                eigen_vectors(3,1:2) = alpha(1:2)
                eigen_vectors(4,1:2) = s(i,1:2)*alpha(1:2)
                eigen_vectors(:,3) = [1.d0,0.d0,-1.d0,0.d0]
                eigen_vectors(:,4) = [1.d0,s(i,4),0.d0,0.d0]
            else
                ! Solve eigenvalue problem
                h_ave(:) = 0.5d0 * (h_l(:) + h_r(:))
                u_ave(:) = 0.5d0 * (u_l(:) + u_r(:))
                A(1,:) = [0.d0,1.d0,0.d0,0.d0]
                A(2,:) = [-u_ave(1)**2 + g*h_ave(1),2*u_ave(1),g*h_ave(1),0.d0]
                A(3,:) = [0.d0,0.d0,0.d0,1.d0]
                A(4,:) = [g*r*h_ave(2),0.d0,-u_ave(2)**2 + g*h_ave(2),2*u_ave(2)]
                call dgeev('N','V',4,A,4,real_evalues,imag_evalues,empty,1,eigen_vectors,4,work,lwork,info)
                if (info < 0) then
                    info = -info
                    print "(a,i1,a)","The ",info,"th argument had an illegal value."
                    stop
                else if (info > 0) then
                    print "(a)","The QR algorithm failed to compute all the"
                    print "(a)","eigenvalues, and no eigenvectors have been"
                    print "(a,i1,a)","computed; elements",i,"+1:4 of WR and WI"
                    print "(a)","contain eigenvalues which have converged."
                    stop
                endif
                do j=1,4
                    eigen_vectors(:,j) = eigen_vectors(:,j) / eigen_vectors(1,j)
                    if (imag_evalues(j) > 0.d0) then
                        print "(a,i1,a,d16.8)","Imaginary eigenvalue(",j,") > 0.0",imag_evalues(j)
                        stop
                    endif
                    s(i,j) = real_evalues(j)
                enddo
!                 if (dry_state_r(2).and.(.not.dry_state_l(2))) then
!                     print *,"Eigenspeeds = "
!                     print "(4d16.8)",(s(i,j),j=1,4)
!                     print *,"Eigenvectors = "
!                     do m=1,4
!                         print "(4d16.8)",(eigen_vectors(m,j),j=1,4)
!                     enddo
!                 endif
            endif
        else
            stop "Invalid eigensystem method requested, method = (1,4)."
        endif
        
        ! ====================================================================
        ! Calculate flux vector to be projected onto e-space
        h_ave(:) = 0.5d0 * (h_r(:) + h_l(:))
        
        ! No dry state
        if ((.not.dry_state_r(2)).and.(.not.dry_state_l(2))) then
            do j=1,2
                flux_r(2*j-1) = rho(j) * hu_r(j)
                flux_l(2*j-1) = rho(j) * hu_l(j)
                flux_r(2*j) = rho(j) * h_r(j) * u_r(j)**2 + 0.5d0 * g * rho(j) * h_r(j)**2
                flux_l(2*j) = rho(j) * h_l(j) * u_l(j)**2 + 0.5d0 * g * rho(j) * h_l(j)**2
            enddo
            flux_r(4) = flux_r(4) + g * h_r(1) * h_r(2) * rho(1)
            flux_l(4) = flux_l(4) + g * h_l(1) * h_l(2) * rho(1)
            
            delta = flux_r - flux_l
        
            ! Note that h_ave include rho values in it
            momentum_transfer = g * rho(1) * h_ave(1) * (h_r(2) - h_l(2))
        
            delta(2) = delta(2) + momentum_transfer + g * rho(1) * h_ave(1) * (b_r - b_l)
            delta(4) = delta(4) - momentum_transfer + g * rho(2) * h_ave(2) * (b_r - b_l)
        ! Right dry state
        else if (dry_state_r(2).and.(.not.dry_state_l(2))) then
            delta(1) = rho(1) * (hu_r(1) - hu_l(1))
            
            eta_r(2) = b_r
            eta_l(2) = h_l(2) + b_l
            
            flux_r(2) = rho(1) * h_r(1) * u_r(1)**2 + 0.5d0 * g * rho(1) * h_r(1)**2
            flux_l(2) = rho(1) * h_l(1) * u_l(1)**2 + 0.5d0 * g * rho(1) * h_l(1)**2
            delta(2) = flux_r(2) - flux_l(2) + g * rho(1) * h_ave(1) * (eta_r(2) - eta_l(2))
            
            h_r(2) = h_l(2)
            hu_r(2) = -hu_l(2)
            u_r(2) = -u_l(2)

            delta(3) = rho(2) * (hu_r(2) - hu_l(2))                                
            flux_r(4) = rho(2) * h_r(2) * u_r(2)**2 + 0.5d0 * g * rho(2) * h_r(2)**2
            flux_l(4) = rho(2) * h_l(2) * u_l(2)**2 + 0.5d0 * g * rho(2) * h_l(2)**2
            delta(4) = flux_r(4) - flux_l(4)
            
        ! Left dry state
        else if (dry_state_l(2).and.(.not.dry_state_r(2))) then
            stop "Not implemented yet..."
        ! Dry state
        else
            delta(1) = rho(1) * (hu_r(1) - hu_l(1))
            flux_r(2) = rho(1) * (h_r(1) * u_r(1)**2 + 0.5d0 * g * h_r(1)**2)
            flux_l(2) = rho(1) * (h_l(1) * u_l(1)**2 + 0.5d0 * g * h_l(1)**2)
            delta(2) = flux_r(2) - flux_l(2) + g * rho(1) * h_ave(1) * (b_r - b_l)
            delta(3:4) = 0.d0
        endif
        
        ! Wind forcing
        wind_speed = 0.5d0 * (w_l + w_r)
        tau = wind_drag(wind_speed) * rho_air * wind_speed
        delta(2) = delta(2) - tau * wind_speed


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
                print *, "Location (i) = (",i,")"
                print *, "Dry states, L=",dry_state_l(2)," R=",dry_state_r(2)
                print *, "h_l(2) = ",h_l(2)," h_r(2) = ",h_r(2)
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
!             print "(a,4d16.8)","Beta = ",(beta(j),j=1,4)
!             print "(a,4d16.8)","  s = ",(s(i,j),j=1,4)
!         endif
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
        do mw=1,mwaves
            if (s(i,mw) > 0.d0) then
                apdq(i,:) = apdq(i,:) + fwave(i,:,mw)
            else                                     
                amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
            endif
        enddo
    enddo

end subroutine rp1