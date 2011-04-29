subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,aux1,aux2,aux3,ilr,asdq,bmasdq,bpasdq)
! ============================================================================
!  Solves transverse Riemann problem for the multilayer shallow water 
!  equations in 2D with topography and wind forcing:
!    (h_1)_t + (h_1 u_1)_x + (h_1 v_1)_y = 0
!    (h_1 u_1)_t + (h_1 u_1^2 + 1/2 g h_1^2)_x + (h_1 u_1 v_1)_y = -gh_1(r(h_2)_x + B_x)
!    (h_1 v_1)_t + (h_1 u_1 v_1)_x + (h_1 v_1^2 + 1/2 g h_1^2)_y = -gh_1(r(h_2)_y + B_y)
!    (h_2)_t + (h_2 u_2)_x + (h_2 v_2)_y = 0
!    (h_2 u_2)_t + (h_2 u_2^2 + 1/2 g h_2^2)_x + (h_2 u_2 v_2)_y = -gh_2(h_1 + B)_x + Tau |W| W_x
!    (h_2 v_2)_t + (h_2 u_2 v_2)_x + (h_2 v_2^2 + 1/2 g h_2^2)_y = -gh_2(h_1 + B)_y + Tau |W| W_y
!
!  On input, ql contains the state vector at the left edge of each cell and qr
!  contains the state vector at the right edge of each cell
!
!           |            |          
!    qr(i-1)|ql(i)  qr(i)|ql(i+1)   
! ----------|------------|----------
!    i-1          i         i+1
!
!  The i-1/2 Riemann problem has left state qr(i-1) and right state ql(i)
!
!  If ixy == 1 then the sweep direction is x, ixy == 2 implies the y direction
!  if ilr == 1 then 
!
!  Kyle T. Mandli (10-13-2010)
! ============================================================================

    use geoclaw_module
    use hurricane_module
    use multilayer_module

    implicit none

    ! Input arguments
    integer, intent(in) :: ixy,maxm,meqn,mwaves,mbc,mx,ilr
    double precision, dimension(1-mbc:maxm+mbc,meqn), intent(in) :: ql,qr
    double precision, dimension(1-mbc:maxm+mbc,meqn), intent(in) :: asdq
    double precision, dimension(1-mbc:maxm+mbc,*), intent(in) :: aux1,aux2,aux3
    
    ! Ouput
    double precision, dimension(1-mbc:maxm+mbc,meqn), intent(out) :: bmasdq,bpasdq
    
    ! Local storage
    integer :: i,j,mw,mu,mu_t,nv,nv_t,info
    double precision :: g,dxdcm,dxdcp,kappa,total_depth,mult_depth
    double precision, dimension(2) :: h,hu,hv,u,v
    double precision, dimension(3) :: b,wx,wy
    double precision, dimension(6) :: s,delta,pivot
    double precision, dimension(6,6) :: eig_vec,A
    
    ! Output array intializations
    bmasdq = 0.d0
    bpasdq = 0.d0

    ! Convenience
    g = grav

    ! Normal and transverse sweep directions
    if (ixy == 1) then
	    mu = 2
	    nv = 3
	    mu_t = 3
	    nv_t = 2
    else
	    mu = 3
	    nv = 2
	    mu_t = 2
	    nv_t = 3
    endif

    ! ========================================================================
    ! Loop over each cell and decompose fluctuations into transverse waves
    ! ========================================================================
    do i=2-mbc,mx+mbc        
        ! Parse states and set appropriate zeros
        ! Note that the "u-direction" is the direction of sweeping which 
        ! could actually be the x or y-directions depending on ixy
        ! Split into either left or right state
        if (ilr == 1) then
            h(1) = qr(i-1,1)
            hu(1) = qr(i-1,mu)
            hv(1) = qr(i-1,nv)
            h(2) = qr(i-1,4)
            hu(2) = qr(i-1,mu+3)
            hv(2) = qr(i-1,nv+3)
            
            ! Auxillary arrays
            b = [aux1(i-1,1), aux2(i-1,1), aux3(i-1,1)]
            wx = [aux1(i-1,mu+2), aux2(i-1,mu+2), aux3(i-1,mu+2)]
            wy = [aux1(i-1,nv+2), aux2(i-1,nv+2), aux3(i-1,nv+2)]
        else
            h(1) = ql(i,1)
            hu(1) = ql(i,mu)
            hv(1) = ql(i,nv)
            h(2) = ql(i,4)
            hu(2) = ql(i,mu+3)
            hv(2) = ql(i,nv+3)
            
            ! Auxillary arrays
            b = [aux1(i,1), aux2(i,1), aux3(i,1)]
            wx = [aux1(i,mu+2), aux2(i,mu+2), aux3(i,mu+2)]
            wy = [aux1(i,nv+2), aux2(i,nv+2), aux3(i,nv+2)]
        endif
        
        ! Apply some dry tolerance limiting and check for bad Riemann problems
        ! and extract velocities
        do j=1,2
            if (h(j) < 0.d0) then
                print "(a,i3,a,2d16.8)","Negative depth value: h(",i,") = ",h(j)
                stop
            endif
            
            if (h(j) < drytolerance) then
                h(j) = 0.d0
                hu(j) = 0.d0
                hv(j) = 0.d0
                u(j) = 0.d0
                v(j) = 0.d0
            else
                u(j) = hu(j) / h(j)
                v(j) = hv(j) / h(j)
            endif
        enddo
        
        ! Convenience variables
        total_depth = sum(h)
        mult_depth = product(h)
        
        ! ====================================================================
        !  Calculate transverse wave speeds
        kappa = (u(1) - u(2))**2 / (g*one_minus_r*total_depth)
        if (kappa > 1.d0) then
            print "(a,2d16.8)","Hyperbolicity may have failed, kappa = ", kappa
!             print "(a,2d16.8)","(dx,dy) = ",dxcom,dycom
            stop
        endif
        ! Approximation based on linearized problem with u_1 = u_2 = v_1 = v_2 = 0
        s(1) = - sqrt(g*total_depth) + 0.5d0 * mult_depth/ total_depth**(3/2) * one_minus_r
        s(2) = - sqrt(g * mult_depth / total_depth * one_minus_r)
        s(3:4) = v
        s(5) = sqrt(g * mult_depth / total_depth * one_minus_r)
        s(6) = sqrt(g*total_depth) - 0.5d0 * mult_depth / total_depth**(3/2) * one_minus_r
            
        
        ! ====================================================================
        !  Calculate eigenvector matrix
        eig_vec(1,:) = [1.d0,1.d0,0.d0,0.d0,1.d0,1.d0]
        eig_vec(mu_t,:) = [s(1),s(2),0.d0,0.d0,s(5),s(6)]
        eig_vec(nv_t,:) = [v(1),v(1),1.d0,0.d0,v(1),v(1)]
        eig_vec(4,1:2) = ((s(1:2)-u(1))**2 - g*h(1)) / (r*g*h(1))
        eig_vec(4,3:4) = 0.d0
        eig_vec(4,5:6) = ((s(5:6)-u(1))**2 - g*h(1)) / (r*g*h(1))
        eig_vec(mu_t+3,:) = s(:) * eig_vec(4,:)
        eig_vec(nv_t+3,1:2) = v(2) * eig_vec(4,1:2)
        eig_vec(nv_t+3,3:4) = [0.d0,1.d0]
        eig_vec(nv_t+3,5:6) = v(2) * eig_vec(4,5:6)
        

        ! ====================================================================
        !  Project asdq onto the eigenspace - Use LAPACK's dgesv routine
        !    N - (int) - Number of linear equations (6)
        !    NRHS - (int) - Number of right hand sides (1)
        !    A - (dp(6,6)) - Coefficient matrix, in this case eig_vec
        !    LDA - (int) - Leading dimension of A (6)
        !    IPIV - (int(N)) - Pivot indices
        !    B - (dp(LDB,NRHS)) - RHS of equations (delta)
        !    LDB - (int) - Leading dimension of B (6)
        !    INFO - (int) - Status of result
        !  Note that the solution (betas) are in delta after the call
        A = eig_vec ! We need to do this as the return matrix is modified and
                    ! we have to use eig_vec again to compute fwaves
        call dgesv(6,1,A,6,pivot,delta,6,info)
        if (.not.(info == 0)) then
            print "(a,i3)","Error solving R beta = delta,",info
            print "(a,i3)","Location: ",i
            print "(a,6d16.8)","Eigenspeeds: ",s(:)
            print "(a)","Eigenvectors:"
            do j=1,6
                print "(6d16.8)",eig_vec(j,:)
            enddo
            stop
        endif

        ! Need to use speeds?
        
        ! ====================================================================
        !  Determine transverse fluctuations
        do mw=1,mwaves
            if (s(mw) > 0.d0) then
                bpasdq(i,:) = bpasdq(i,:) + dxdcp * s(mw) * delta(mw) * eig_vec(:,mw)
            else
                bmasdq(i,:) = bmasdq(i,:) + dxdcm * s(mw) * delta(mw) * eig_vec(:,mw)
            endif
        enddo
    enddo
    
end subroutine rpt2