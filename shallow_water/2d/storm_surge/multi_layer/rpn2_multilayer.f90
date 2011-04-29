subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! ============================================================================
!  Solves normal Riemann problem for the multilayer shallow water equations in
!  2D with topography and wind forcing:
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
!
!  Kyle T. Mandli (10-11-2010)
! ============================================================================

    use geoclaw_module
    use hurricane_module
    use multilayer_module

    implicit none

    ! Input arguments
    integer, intent(in) :: ixy,maxm,meqn,mwaves,mbc,mx
    double precision, dimension(1-mbc:maxm+mbc,meqn), intent(in) :: ql,qr
    double precision, dimension(1-mbc:maxm+mbc,*), intent(in) :: auxl,auxr

    ! Output arguments
    double precision, dimension(1-mbc:maxm+mbc, meqn, mwaves), intent(out) :: fwave
    double precision, dimension(1-mbc:maxm+mbc, mwaves), intent(out) :: s
    double precision, dimension(1-mbc:maxm+mbc, meqn), intent(out) :: apdq, amdq

    ! Local variables
    integer :: i,j,m,mw,mu,nv,k,maxiter,info
    double precision :: g,wind_speed,tau,dxdc
    double precision, dimension(2) :: h_l,h_r,hu_l,hu_r,hv_l,hv_r
    double precision, dimension(2) :: u_l,u_r,v_l,v_r,advected_speed
    double precision :: b_l,b_r,wx,wy,kappa_l,kappa_r
    double precision :: total_depth_l,total_depth_r,mult_depth_l,mult_depth_r
    double precision, dimension(6) :: delta,flux_r,flux_l,pivot
    double precision, dimension(6,6) :: eig_vec,A

    ! Common block variables
    integer :: mcapa,icom,jcom
    double precision :: dtcom,dxcom,dycom,tcom

    common /cmcapa/ mcapa
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

    external dgesv

    ! Initialize output variables
    amdq = 0.d0
    apdq = 0.d0
    
    ! Set normal direction
    if (ixy == 1) then
        mu = 2
        nv = 3
    else
        mu = 3
        nv = 2
    endif

    g = grav
    
    ! ========================================================================
    ! Loop through Riemann problems
    ! ========================================================================
    do i=2-mbc,mx+mbc
        ! Parse states and set appropriate zeros
        ! Note that the "u-direction" is the direction of sweeping which 
        ! could actually be the x or y-directions depending on ixy
        h_l(1) = qr(i-1,1)
        hu_l(1) = qr(i-1,mu)
        hv_l(1) = qr(i-1,nv)
        
        h_l(2) = qr(i-1,4)
        hu_l(2) = qr(i-1,mu+3)
        hv_l(2) = qr(i-1,nv+3)
        
        h_r(1) = ql(i,1)
        hu_r(1) = ql(i,mu)
        hv_r(1) = ql(i,nv)

        h_r(2) = ql(i,4)
        hu_r(2) = ql(i,mu+3)
        hv_r(2) = ql(i,nv+3)
        
        b_l = auxr(i-1,1)
        b_r = auxl(i,1)
        wx = 0.5d0 * (auxr(i-1,mu+2) + auxl(i,mu+2))
        wy = 0.5d0 * (auxr(i-1,nv+2) + auxl(i,nv+2))

        ! Apply some dry tolerance limiting and check for bad Riemann problems
        ! and extract velocities
        do j=1,2
            if (h_l(j) < 0.d0) then
                print "(a,i3,a,2d16.8)","Negative depth value: hl(",i,") = ",h_l
                stop
            endif
            if (h_r(j) < 0.d0) then
                print "(a,i3,a,2d16.8)","Negative depth value: hr(",i,") = ",h_r
                stop
            endif
            
            if (h_l(j) < drytolerance) then
                h_l(j) = 0.d0
                hu_l(j) = 0.d0
                hv_l(j) = 0.d0
                u_l(j) = 0.d0
                v_l(j) = 0.d0
            else
                u_l(j) = hu_l(j) / h_l(j)
                v_l(j) = hv_l(j) / h_l(j)
            endif
            if (h_r(j) < drytolerance) then
                h_r(j) = 0.d0
                hu_r(j) = 0.d0
                hv_r(j) = 0.d0
                u_r(j) = 0.d0
                v_r(j) = 0.d0
            else
                u_r(j) = hu_r(j) / h_r(j)
                v_r(j) = hv_r(j) / h_r(j)
            endif
        enddo
        
        ! Convenience variables
        total_depth_l = sum(h_l)
        total_depth_r = sum(h_r)
        mult_depth_l = product(h_l)
        mult_depth_r = product(h_r)
        
        ! Figure out which speed to use for transverse direction
        do j=1,2
            if ((u_l(j) < 0.d0).and.(u_r(j) < 0.d0)) then
                advected_speed(j) = u_l(j)
            else if ((u_l(j) > 0.d0).and.(u_r(j) > 0.d0)) then
                advected_speed(j) = u_r(j)
            else
!                 print *,"Transonic Riemann problem detected!"
!                 print "(a,d16.8,a,d16.8,a)","(u_l,u_r) = (",u_l(j),",",u_r(j),")"
                advected_speed(j) = 0.5d0 * (u_l(j) + u_r(j))
            endif
        enddo
        
        ! ====================================================================
        !  Calculate wave speeds
        ! ====================================================================
        kappa_l = (u_l(1) - u_l(2))**2 / (g*one_minus_r*total_depth_l)
        kappa_r = (u_r(1) - u_r(2))**2 / (g*one_minus_r*total_depth_r)
        if ((kappa_l > 1.d0).or.(kappa_r) > 1.d0) then
            print "(a,2d16.8)","Hyperbolicity may have failed, kappa = ", kappa_l,kappa_r
            print "(a,2d16.8)","(dx,dy) = ",dxcom,dycom
            stop
        endif
        ! Approximation based on linearized problem with u_1 = u_2 = v_1 = v_2 = 0
        s(i,1) = - sqrt(g*total_depth_l) + 0.5d0 * mult_depth_l/ total_depth_l**(3/2) * one_minus_r
        s(i,2) = - sqrt(g * mult_depth_l / total_depth_l * one_minus_r)
        s(i,3:4) = advected_speed
        s(i,5) = sqrt(g * mult_depth_r / total_depth_r * one_minus_r)
        s(i,6) = sqrt(g*total_depth_r) - 0.5d0 * mult_depth_r / total_depth_r**(3/2) * one_minus_r
            
        ! ====================================================================
        ! Calculate eigenvector matrix
        eig_vec(1,:) = [1.d0,1.d0,0.d0,0.d0,1.d0,1.d0]
        eig_vec(mu,:) = [s(i,1),s(i,2),0.d0,0.d0,s(i,5),s(i,6)]
        eig_vec(nv,:) = [v_l(1),v_l(1),1.d0,0.d0,v_r(1),v_r(1)]
        eig_vec(4,1:2) = ((s(i,1:2)-u_l(1))**2 - g*h_l(1)) / (r*g*h_l(1))
        eig_vec(4,3:4) = 0.d0
        eig_vec(4,5:6) = ((s(i,5:6)-u_r(1))**2 - g*h_r(1)) / (r*g*h_r(1))
        eig_vec(mu+3,:) = s(i,:) * eig_vec(4,:)
        eig_vec(nv+3,1:2) = v_l(2) * eig_vec(4,1:2)
        eig_vec(nv+3,3:4) = [0.d0,1.d0]
        eig_vec(nv+3,5:6) = v_r(2) * eig_vec(4,5:6)

        ! ====================================================================
        !  Compute jump in fluxes
        flux_r(1) = hu_r(1)
        flux_r(mu) = hu_r(1)**2 / h_r(1) + 0.5d0 * g * h_r(1)**2
        flux_r(nv) = hu_r(1) * hv_r(1) / h_r(1)
        flux_r(4) = hu_r(2)
        flux_r(mu+3) = hu_r(2)**2 / h_r(2) + 0.5d0 * g * h_r(2)**2
        flux_r(nv+3) = hu_r(2) * hv_r(2) / h_r(2)
        
        flux_l(1) = hu_l(1)
        flux_l(mu) = hu_l(1)**2 / h_l(1) + 0.5d0 * g * h_l(1)**2
        flux_l(nv) = hu_l(1) * hv_l(1) / h_l(1)
        flux_l(4) = hu_l(2)
        flux_l(mu+3) = hu_l(2)**2 / h_l(2) + 0.5d0 * g * h_l(2)**2
        flux_l(nv+3) = hu_l(2) * hv_l(2) / h_l(2)
        
        delta = flux_r - flux_l
        
        ! Add in non-conservative product and bathymetry and wind source terms
        delta(mu) = delta(mu) + 0.5d0 * g * (h_l(1) + h_r(1)) * (r * (h_r(2)-h_l(2)) + b_r - b_l)
        delta(mu+3) = delta(mu+3) + 0.5d0 * g * (h_l(2) + h_r(2)) * (h_r(1) - h_l(1) + b_r - b_l)
!         wind_speed = sqrt(wx**2 + wy**2)
!         tau = wind_drag(wind_speed) * rho_air * wind_speed
!         delta(mu+3) = delta(mu+3) - tau * wx
        
        ! ====================================================================
        !  Project jump in fluxes - Use LAPACK's dgesv routine
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
            print "(a,6d16.8)","Eigenspeeds: ",s(i,:)
            print "(a)","Eigenvectors:"
            do j=1,6
                print "(6d16.8)",eig_vec(j,:)
            enddo
            stop
        endif

        ! ====================================================================
        ! Compute fwaves
        forall(mw=1:mwaves)
            fwave(i,:,mw) = eig_vec(:,mw) * delta(mw)
        endforall
        
    enddo
    ! ========================================================================
    
    ! ====================================================================
    ! Capacity for mapping from latitude longitude to physical space
    if (mcapa > 0) then
        do i=2-mbc,mx+mbc
            if (ixy == 1) then
                dxdc=(Rearth*pi/180.d0)
            else
                dxdc=auxl(i,3)
            endif

            do mw=1,mwaves
    	        s(i,mw)=dxdc*s(i,mw)
    	        fwave(i,:,mw)=dxdc*fwave(i,:,mw)
            enddo
        enddo
    endif

    ! ========================================================================
    !  Compute fluctuations 
    do i=2-mbc,mx+mbc
        do mw=1,mwaves
            if (s(i,mw) > 0.d0) then
                apdq(i,:) = apdq(i,:) + fwave(i,:,mw)
            else
                amdq(i,:) = amdq(i,:) + fwave(i,:,mw)
            endif
        enddo
    enddo

end subroutine rpn2
