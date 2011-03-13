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
    integer :: i,j,m,mw
    double precision :: h_l(2),hu_l(2),h_r(2),hu_r(2),u_l(2),u_r(2),b_l,b_r
    double precision :: h_star(2),temp_s(4),temp_fwave(4,4)

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
        ! Extract state variables, here h and hu have values corresponding to
        ! which layer they are in
        h_l(1) = qr(i-1,1)
        hu_l(1) = qr(i-1,2)
        h_r(1) = ql(i,1)
        hu_r(1) = ql(i,2)
        
        h_l(2) = qr(i-1,3)
        hu_l(2) = qr(i-1,4)
        u_l(2) = hu_l(2) / h_l(2) ! Assume top layer is not zero
        h_r(2) = ql(i,3)
        hu_r(2) = ql(i,4)
        u_l(2) = hu_r(2) / h_r(2) ! Assume top layer is not zero
        
        b_l = auxr(i-1,1)
        b_r = auxl(i,1)
        
        ! Assuming right now that top layer does not goto 0
        if ((h_l(2) < DRY_TOLERANCE).or.(h_r(2) < DRY_TOLERANCE)) then
            print "(a,i4)", "Dry state in top layer found at i=",i,", aborting."
            stop
        endif
        
        ! Check for bottom layer dry state
        if ((h_l(1) < DRY_TOLERANCE).or.(h_r(1) < DRY_TOLERANCE)) then
            ! Left side dry, solve wall Riemann problem for right going waves
            if ((h_l(1) < DRY_TOLERANCE).and.(h_r(1) > DRY_TOLERANCE)) then
                stop "Not implemented yet!"
            ! Rigth side dry, solve wall Riemann problem for left going waves
            else if ((h_r(1) < DRY_TOLERANCE).and.(h_l(1) > DRY_TOLERANCE)) then
                ! Solve as bottom layer wall problem
                u_l(1) = hu_l(1) / h_l(1)
                
                h_r(1) = h_l(1)
                hu_r(1) = -hu_l(1)
                u_r(1) = u_l(1)
                
                ! Solve multilayer equations
                call rp_multilayer(i,h_l,h_r,hu_l,hu_r,u_l,u_r,b_l,b_r,auxl(i,2),auxr(i,2),h_star,temp_s,temp_fwave)
                s(i,:) = temp_s
                fwave(i,:,:) = temp_fwave
                
                ! Throw out 3rd wave
                fwave(i,:,3) = 0.d0
                s(i,3) = 0.d0
                
                ! Check to see if h_l > jump in B
                if (h_star(1) > abs(b_r - b_l)) then
                    ! Dry cell has become wet, must recalculate
                    stop "Not quite yet..."
                endif
                
            ! Both sides are dry, do nothing and solve swe problem
            else
                h_l(1) = 0.d0
                hu_l(1) = 0.d0
                h_r(1) = 0.d0
                hu_r(1) = 0.d0
                
                ! Solution
                s(i,2:3) = 0.d0  ! Internal waves are zero
                fwave(i,3:4,2:3) = 0.d0  ! External waves
                fwave(i,1:2,:) = 0.d0
                
                ! Solve single layer swe equations
                call rp_singlelayer(i,h_l(2),h_r(2),hu_l(2),hu_r(2),b_l,b_r,h_star(1),temp_s,temp_fwave)
                s(i,:) = temp_s
                fwave(i,:,:) = temp_fwave
            endif
        else
            u_l(1) = hu_l(1) / h_l(1)
            u_r(1) = hu_r(1) / h_r(1)
            
            ! Solve multilayer equations
            call rp_multilayer(i,h_l,h_r,hu_l,hu_r,u_l,u_r,b_l,b_r,auxl(i,2),auxr(i,2),h_star,temp_s,temp_fwave)
    
            if ((i>48).and.(i<52)) then
                print "(a)","====== rp1 1 ================"
                do mw=1,4
                    print "('fwave ',i1,'=',4d16.8)",mw,temp_fwave(:,mw)
                enddo
                print "(a)","============================="
            endif
            
            s(i,:) = temp_s
            fwave(i,:,:) = temp_fwave
            
            if ((i>48).and.(i<52)) then
                print "(a)","====== rp1 2 ================"
                do mw=1,4
                    print "('fwave ',i1,'=',4d16.8)",mw,fwave(i,1:4,mw)
                enddo
                print "(a)","============================="
            endif
        endif
        
    enddo
    
    ! Calculate amdq and apdq
    do i=2-mbc,mx+mbc
        do mw=1,mwaves
            if (s(i,mw) > 0.d0) then
                apdq(i,:) = apdq(i,:) + fwave(i,:,mw)
            else
                amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
            endif
        enddo
    enddo
!     forall(i=2-mbc:mx+mbc,m=1:meqn,s(i,mw) > 0.d0)
!         apdq(i,m) = apdq(i,m) + sum(fwave(i,m,:))
!     end forall
!     forall(i=2-mbc:mx+mbc,m=1:meqn,s(i,mw) <= 0.d0)
!         amdq(i,m) = amdq(i,m) + sum(fwave(i,m,:))
!     end forall
!     do mw=1,mwaves
!         if (s(i,mw) > 0.d0) then
!             apdq(i,:) = apdq(i,:) + fwave(i,:,mw)
!         else
!             amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
!         endif
!     enddo

end subroutine rp1

subroutine rp_singlelayer(i,h_l,h_r,hu_l,hu_r,b_l,b_r,h_star,s,fwave)

    use parameters_module

    implicit none
    integer, intent(in) :: i
    double precision, intent(in) :: h_l,h_r,hu_l,hu_r,b_l,b_r
    double precision, intent(inout) :: h_star,s(4),fwave(4,4)
    
    integer :: mw
    double precision :: h_bar,u_bar
    double precision, dimension(2) :: delta,alpha
    
    s = 0.d0
    fwave = 0.d0
    
    h_bar = 0.5d0 * (h_l + h_r)
    u_bar = (hu_l / sqrt(h_l) + hu_r / sqrt(h_r)) / (sqrt(h_l) + sqrt(h_r))
    
    s(1) = u_bar - sqrt(g * h_bar)
    s(4) = u_bar + sqrt(g * h_bar)
    
    delta(1) = hu_r - hu_l
    delta(2) = hu_r**2 / h_r + 0.5d0*g*h_r**2 - (hu_l**2 / h_l + 0.5d0*g*h_l**2)
    delta(2) = delta(2)  + 0.5d0 * g * (h_l + h_r) * (b_r - b_l)
    
    alpha(1) = (delta(2) - s(2) * delta(1)) / (s(1) - s(2))
    alpha(2) = delta(1) - alpha(1)

    fwave(:,1) = alpha(1) * [0.d0,0.d0,1.d0,s(1)]
    fwave(:,4) = alpha(2) * [0.d0,0.d0,1.d0,s(4)]
    
end subroutine rp_singlelayer

subroutine rp_multilayer(i,h_l,h_r,hu_l,hu_r,u_l,u_r,b_l,b_r,w_l,w_r,h_star,s,fwave)

    use parameters_module

    implicit none
    
    ! State variables
    integer, intent(in) :: i
    double precision, intent(in), dimension(2) :: h_l,h_r,hu_l,hu_r,u_l,u_r
    double precision, intent(in) :: b_l,b_r,w_l,w_r
    
    ! Output
    double precision, intent(inout) :: h_star(2),s(4), fwave(4,4)
    
    ! Local variables
    integer :: j,mw,info,ipiv(4)
    double precision :: total_depth_l,total_depth_r,mult_depth_l,mult_depth_r
    double precision :: e2(4),s2(4,2),kappa_l,kappa_r
    double precision :: flux_l(4),flux_r(4),delta(4),eig_vec(4,4),A(4,4)
    double precision :: wind_speed,tau
    
    ! Convenience variables
    total_depth_l = sum(h_l)
    total_depth_r = sum(h_r)
    mult_depth_l = product(h_l)
    mult_depth_r = product(h_r)
    
    ! ============================================================================
    !  Calculate wave speeds
    ! ============================================================================
    ! Check hyperbolicity condition
    kappa_l = (u_l(1)-u_l(2))**2 / (g*one_minus_r*total_depth_l)
    kappa_r = (u_r(1)-u_r(2))**2 / (g*one_minus_r*total_depth_r)
    if ((kappa_l > 1.d0).or.(kappa_r > 1.d0)) then
        print "(a,2d16.8)","Hyperbolicity may have failed, kappa = ",kappa_l,kappa_r
    endif
    ! Approximation based on (u_1 - u_2) being small, commonly found in
    ! literature
    if (eigenvalue_type == 1) then
        s(1) = (h_l(1)*u_l(1)+h_l(2)*u_l(2)) / total_depth_l - sqrt(g*total_depth_l)
        s(2) = (h_l(2)*u_l(1)+h_l(1)*u_l(2)) / total_depth_l - sqrt(g*one_minus_r*mult_depth_l/total_depth_l * (1-(u_l(1)-u_l(2))**2/(g*one_minus_r*total_depth_l)))
        s(3) = (h_r(2)*u_r(1)+h_r(1)*u_r(2)) / total_depth_r + sqrt(g*one_minus_r*mult_depth_r/total_depth_r * (1-(u_r(1)-u_r(2))**2/(g*one_minus_r*total_depth_r)))
        s(4) = (h_r(1)*u_r(1)+h_r(2)*u_r(2)) / total_depth_r + sqrt(g*total_depth_r)
    else
        s2(1,1) = (h_l(1)*u_l(1)+h_l(2)*u_l(2)) / total_depth_l - sqrt(g*total_depth_l)
        s2(2,1) = (h_l(2)*u_l(1)+h_l(1)*u_l(2)) / total_depth_l - sqrt(g*one_minus_r*mult_depth_l/total_depth_l * (1-(u_l(1)-u_l(2))**2/(g*one_minus_r*total_depth_l)))
        s2(3,1) = (h_r(2)*u_r(1)+h_r(1)*u_r(2)) / total_depth_r + sqrt(g*one_minus_r*mult_depth_r/total_depth_r * (1-(u_r(1)-u_r(2))**2/(g*one_minus_r*total_depth_r)))
        s2(4,1) = (h_r(1)*u_r(1)+h_r(2)*u_r(2)) / total_depth_r + sqrt(g*total_depth_r)
    endif

    ! Approximation based on linearized system about u_1 = u_2 = 0
    if (eigenvalue_type == 2) then
        s(1) = - sqrt(g*total_depth_l) + 0.5d0 * mult_depth_l/ total_depth_l**(3/2) * one_minus_r
        s(2) = - sqrt(g * mult_depth_l / total_depth_l * one_minus_r)
        s(3) = sqrt(g * mult_depth_r / total_depth_r * one_minus_r)
        s(4) = sqrt(g*total_depth_r) - 0.5d0 * mult_depth_r / total_depth_r**(3/2) * one_minus_r
    else
        s2(1,1) = - sqrt(g*total_depth_l) + 0.5d0 * mult_depth_l/ total_depth_l**(3/2) * one_minus_r
        s2(2,1) = - sqrt(g * mult_depth_l / total_depth_l * one_minus_r)
        s2(3,1) = sqrt(g * mult_depth_r / total_depth_r * one_minus_r)
        s2(4,1) = sqrt(g*total_depth_r) - 0.5d0 * mult_depth_r / total_depth_r**(3/2) * one_minus_r
    endif
    
    ! Approximation as above except eta_1 = eta_2 = 0 as well, leads to a
    ! time independent system
!         s(i,1) = - sqrt(g*abs(b_l)) + abs(H) * (abs(b_l)-abs(H)) / (2.d0 * abs(b_l)**(3/2))*one_minus_r
!         s(i,2) = - sqrt(g*abs(H)*(abs(b_l)-abs(H)) / abs(b_l)*one_minus_r)
!         s(i,3) = sqrt(g*abs(H)*(abs(b_r)-abs(H)) / abs(b_r) * one_minus_r)
!         s(i,4) = sqrt(g*abs(b_r)) - abs(H) * (abs(b_r)-abs(H)) / (2.d0 * abs(b_r)**(3/2))*one_minus_r
!         s2(1,2) = - sqrt(g*abs(b_l)) + abs(H) * (abs(b_l)-abs(H)) / (2.d0 * abs(b_l)**(3/2))*one_minus_r
!         s2(2,2) = - sqrt(g*abs(H)*(abs(b_l)-abs(H)) / abs(b_l)*one_minus_r)
!         s2(3,2) = sqrt(g*abs(H)*(abs(b_r)-abs(H)) / abs(b_r) * one_minus_r)
!         s2(4,2) = sqrt(g*abs(b_r)) - abs(H) * (abs(b_r)-abs(H)) / (2.d0 * abs(b_r)**(3/2))*one_minus_r  

    ! Check differences
!         print "(a,i3,a)","== i = ",i," =============================="
!         print "(a,i3,a,4d16.8)"," s(",i,",:) = ",s(i,:)
!         print "(a,4d16.8)"," s2(:,1) = ",s2(:,1)
!         print "(a,4d16.8)"," s2(:,2) = ",s2(:,2)
!         print "(a)","========================================="

    ! ============================================================================
    !  Calculate eigenvector matrix
    ! ============================================================================
    ! Calculate eigenvectors corresponding to left going waves based on 
    ! left state and similarly for right going waves
    eig_vec(1,:) = 1.d0
    eig_vec(2,:) = s(:)
    
    ! Need to use coefficient that includes only top layer
!         e2(1:2) = (s(i,1:2)**2 - g*h_l(1)) / (g*r*h_l(1))
!         e2(3:4) = (s(i,3:4)**2 - g*h_r(1)) / (g*r*h_r(1))
    eig_vec(3,1:2) = g * h_l(2) / (s(1:2)**2 - g*h_l(2))
!         print "(a,i4,a)","====== i=",i," ======"
!         print *,"g h_2 = ",g*h_r(2)
!         print *,"s(3,4) = ",(s(i,m),m=3,4)
!         print *,"s(3,4)**2 - g h2 = ",(s(i,m)**2 - g*h_r(2),m=3,4)
!         print *,g * h_r(2) / (s(i,3:4)**2 - g*h_r(2))
    eig_vec(3,3:4) = g * h_r(2) / (s(3:4)**2 - g*h_r(2))
    
    ! Second form, deals with only top layer
!         eig_vec(3,1:2) = g * h_l(2) / (s(i,1:2)**2 - g*h_l(2))
!         eig_vec(3,3:4) = g * h_r(2) / (s(i,3:4)**2 - g*h_r(2))
!         e2(1:2) = (s(i,1:2)**2 - g*h_l(1)) / (g*r*h_l(1))
!         e2(3:4) = (s(i,3:4)**2 - g*h_r(1)) / (g*r*h_r(1))
    
    ! Same for both
    eig_vec(4,:) = s(:) * eig_vec(3,:)

    ! Check differences
!         if (abs(sum(e2(:) - eig_vec(3,:))) > 1.d-8) then
!             print *,"Third component disagrees"
!             print *," i=",i
!             print "(4d16.8)",eig_vec(3,:)
!             print "(4d16.8)",e2
!             print "(4d16.8)",e2(:)-eig_vec(3,:)
!             print "(d16.8)",abs(sum(e2(:)-eig_vec(3,:)))
!             stop
!         endif
    
    ! Jump in fluxes
    flux_r(1) = hu_r(1)
    flux_r(2) = hu_r(1)**2 / h_r(1) + 0.5d0 * g * h_r(1)**2
    flux_r(3) = hu_r(2)
    flux_r(4) = hu_r(2)**2 / h_r(2) + 0.5d0 * g * h_r(2)**2
    
    flux_l(1) = hu_l(1)
    flux_l(2) = hu_l(1)**2 / h_l(1) + 0.5d0 * g * h_l(1)**2
    flux_l(3) = hu_l(2)
    flux_l(4) = hu_l(2)**2 / h_l(2) + 0.5d0 * g * h_l(2)**2

    delta = flux_r - flux_l
    
    ! Add in non-conservative product and bathymetry source terms
    delta(2) = delta(2) + 0.5d0 * g * (h_l(1)+h_r(1)) * (r * (h_r(2)-h_l(2)) + b_r - b_l)
    delta(4) = delta(4) + 0.5d0 * g * (h_l(2)+h_r(2)) * (h_r(1)-h_l(1) + b_r - b_l)
    
    ! Wind forcing, acts only on top layer
    wind_speed = 0.5d0 * (w_l + w_r)
    tau = wind_drag(wind_speed) * rho_air * wind_speed
    delta(4) = delta(4) - tau * wind_speed
!     
!     if ((i>48).or.(i<52)) then
!         print "('delta =',4d16.8)",delta
!         print "(a)","============================="
!     endif
    
    ! Wave strengths, stored in delta after this call
    ! call dgesv(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
    ! N = number of linear equations
    ! NRHS = Number of right hand sides
    ! A = Coefficient matrix A
    ! LDA = Leading dimension of A
    ! IPIV = Pivot indices (dimension(N))
    ! B = RHS of equations dimension(LDB,NRHS)
    ! LDB = Leading dimension of B
    ! INFO = status integer
    A = eig_vec
    call dgesv(4,1,A,4,ipiv,delta,4,info)
    if (.not.(info == 0)) then
        print *, "Error solving R beta = delta, ",info
        print *, "Eigen-speeds:",s(:)
        print *, "Eigen-vectors:"
        do j=1,4
            print "(4d16.8)",eig_vec(j,:)
        enddo
        stop
    endif
    
!     if ((i>48).or.(i<52)) then
!         print "('beta = ',4d16.8)", delta
!         print "('eig_vec = ',4d16.8)", eig_vec
!         print "(a)","============================="
!     endif
    
    ! Calculate waves, note that beta is stored in delta
    forall(mw=1:4)
        fwave(:,mw) = eig_vec(:,mw) * delta(mw)
    end forall
    
    if ((i>48).and.(i<52)) then
        print "(a)","====== rp_multilayer ========"
        do mw=1,4
            print "('fwave ',i1,'=',4d16.8)",mw,fwave(:,mw)
        enddo
        print "(a)","============================="
    endif
    
    ! Calculate h_star, middle state
    if ((abs(s(1)) > 0.d0).and.(abs(s(2)) > 0.d0)) then
        h_star(1) = h_l(1) + 1.d0 / s(1) * fwave(1,1) + 1.d0 / s(2) * fwave(1,2)
        h_star(2) = h_l(2) + 1.d0 / s(1) * fwave(3,1) + 1.d0 / s(2) * fwave(3,2)
    else
        h_star = 1d-99
    endif
end subroutine rp_multilayer

