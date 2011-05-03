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

    ! Counters
    integer :: i,j,m,mw,k,maxiter,info
    integer :: n_index,t_index,layer_index
    
    ! Physics
    double precision :: g,dxdc
    
    ! State variables
    double precision, dimension(2) :: h_l,h_r,hu_l,hu_r,hv_l,hv_r
    double precision, dimension(2) :: u_l,u_r,v_l,v_r,advected_speed
    double precision :: b_l,b_r,w_normal,w_transverse,kappa_l,kappa_r
    double precision :: eta_l(2),eta_r(2),h_ave(2),momentum_transfer(2)
    double precision :: h_hat_l(2),h_hat_r(2),gamma_l,gamma_r
    double precision :: flux_transfer_l,flux_transfer_r

    ! Solver variables
    double precision, dimension(6) :: delta,flux_r,flux_l,pivot
    double precision, dimension(6,6) :: eig_vec,A
    double precision :: beta(6),alpha(4)
    logical :: dry_state_l(2), dry_state_r(2)
    
    ! Single layer locals
    integer, parameter :: max_iterations = 1
    double precision :: wall(3),fw(3,3),sw(3),phi_r(2),phi_l(2)
    double precision :: s_l,s_r,s_roe(2),s_E(2),u_hat,c_hat,sm(2)
    double precision :: h_star,h_star_test,h_star_HLL,s_l_test,s_r_test
    logical :: rare(2)

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
        n_index = 2
        t_index = 3
    else
        n_index = 3
        t_index = 2
    endif

    g = grav
    
    ! ========================================================================
    ! Loop through Riemann problems
    ! ========================================================================
    do i=2-mbc,mx+mbc
        dry_state_l = .false.
        dry_state_r = .false.
        
        ! Parse states and set appropriate zeros
        ! Note that the "u-direction" is the direction of sweeping which 
        ! could actually be the x or y-directions depending on ixy
        
        do j=1,2
            layer_index = 3*(j-1)
            h_l(j) = qr(i-1,layer_index+1) / rho(j)
            hu_l(j) = qr(i-1,layer_index+n_index) / rho(j)
            hv_l(j) = qr(i-1,layer_index+t_index) / rho(j)
            
            h_r(j) = ql(i,layer_index+1) / rho(j)
            hu_r(j) = ql(i,layer_index+n_index) / rho(j)
            hv_r(j) = ql(i,layer_index+t_index) / rho(j)
            
            h_hat_l(j) = auxr(i-1,j+6)
            h_hat_r(j) = auxl(i,j+6)
            
            h_ave(:) = 0.5d0 * (h_l(:) + h_r(:))
            
            ! Check for dry states
            if (h_l(j) < drytolerance) then
!                 print *,"==================="
!                 print "(a)","Left side dry"
!                 print "(a,2d16.8)","h_l",h_l
!                 print "(a,2d16.8)","h_r",h_r
                
                dry_state_l(j) = .true.
                hu_l(j) = 0.d0
                hv_l(j) = 0.d0
                u_l(j) = 0.d0
                v_l(j) = 0.d0
            else
                u_l(j) = hu_l(j) / h_l(j)
                v_l(j) = hv_l(j) / h_l(j)
            endif
            if (h_r(j) < drytolerance) then
!                 print *,"==================="
!                 print "(a)","Right side dry"
!                 print "(a,2d16.8)","h_l",h_l
!                 print "(a,2d16.8)","h_r",h_r
                
                dry_state_r(j) = .true.
                hu_r(j) = 0.d0
                hv_r(j) = 0.d0
                u_r(j) = 0.d0
                v_r(j) = 0.d0
            else
                u_r(j) = hu_r(j) / h_r(j)
                v_r(j) = hv_r(j) / h_r(j)
            endif
        enddo
        
        b_l = auxr(i-1,1)
        b_r = auxl(i,1)
        
        ! Calculate wind stress
!         w_normal = 0.5d0 * (auxr(i-1,n_index+2) + auxl(i,n_index+2))
!         w_transverse = 0.5d0 * (auxr(i-1,t_index+2) + auxl(i,t_index+2))
!         wind_speed = sqrt(w_normal**2 + w_transverse**2)
!         tau = wind_drag(wind_speed) * rho_air * wind_speed
!         if (ixy == 1) then
!             wind_stress = dxcom * tau * w_normal
!         else if (ixy == 2) then
!             wind_stress = dycom * tau * w_normal
!         endif

        ! ====================================================================
        ! Dry state Handling
        ! ====================================================================
        ! Single layer case
        if (dry_state_l(2).and.dry_state_r(2)) then
            wall = 1.d0
            
            ! Completely dry cell
            if (dry_state_l(1).and.dry_state_r(1)) then
                s(i,:) = 0.d0
                fwave(i,:,:) = 0.d0
                cycle
            endif
            
            ! Calculate momentum fluxes
            phi_l(1) = 0.5d0 * g * h_l(1)**2 + h_l(1) * u_l(1)**2
            phi_r(1) = 0.5d0 * g * h_r(1)**2 + h_r(1) * u_r(1)**2
             
            ! Check for dry state to right
            if (h_r(1) <= drytolerance) then
                call riemanntype(h_l(1),h_l(1),u_l(1),-u_l(1),h_star, &
                                 sm(1),sm(2),rare(1),rare(2),1,drytolerance,g)
                h_star_test = max(h_l(1),h_star)
                ! Right state should become ghost values that mirror left for wall problem
                if (h_star_test + b_l < b_r) then 
                    wall(2:3)=0.d0
                    h_r(1) = h_l(1)
                    hu_r(1) = -hu_l(1)
                    b_r = b_l
                    phi_r(1) = phi_l(1)
                    u_r(1) = -u_l(1)
                    v_r(1) = v_l(1)
                elseif (h_l(1) + b_l < b_r) then
                    b_r = h_l(1) + b_l
                endif
            ! Check for drystate to left, i.e right surface is lower than left topo
            else if (h_l(1) <= drytolerance) then 
                call riemanntype(h_r(1),h_r(1),-u_r(1),u_r(1),h_star, &
                                 sm(1),sm(2),rare(1),rare(2),1,drytolerance,g)
                h_star_test = max(h_r(1),h_star)
                ! Left state should become ghost values that mirror right
                if (h_star_test + b_r < b_l) then  
                   wall(1:2) = 0.d0
                   h_l(1) = h_r(1)
                   hu_l(1) = -hu_r(1)
                   b_l = b_r
                   phi_l(1) = phi_r(1)
                   u_l(1) = -u_r(1)
                   v_l(1) = v_r(1)
                elseif (h_r(1) + b_r < b_l) then
                   b_l = h_r(1) + b_r
                endif
             endif

             ! Determine wave speeds
             s_l = u_l(1) - sqrt(g*h_l(1)) ! 1 wave speed of left state
             s_r = u_r(1) + sqrt(g*h_r(1)) ! 2 wave speed of right state
             
             u_hat = (sqrt(g*h_l(1))*u_l(1) + sqrt(g*h_r(1))*u_r(1)) &
                        / (sqrt(g*h_r(1))+sqrt(g*h_l(1))) ! Roe average
             c_hat = sqrt(g*0.5d0*(h_r(1)+h_l(1))) ! Roe average
             s_roe(1) = u_hat - c_hat ! Roe wave speed 1 wave
             s_roe(2) = u_hat + c_hat ! Roe wave speed 2 wave

             s_E(1) = min(s_l,s_roe(1)) ! Eindfeldt speed 1 wave
             s_E(2) = max(s_r,s_roe(2)) ! Eindfeldt speed 2 wave

             ! Solve Riemann problem
             call riemann_aug_JCP(max_iterations,3,3,h_l(1),h_r(1),hu_l(1), &
                                    hu_r(1),hv_l(1),hv_r(1),b_l,b_r,u_l(1), &
                                    u_r(1),v_l(1),v_r(1),phi_l(1),phi_r(1), &
                                    s_E(1),s_E(2),drytolerance,g,sw,fw)
            
!             call riemann_fwave(max_iterations,3,3,h_l(1),h_r(1),hu_l(1), &
!                                     hu_r(1),hv_l(1),hv_r(1),b_l,b_r,u_l(1), &
!                                     u_r(1),v_l(1),v_r(1),phi_l(1),phi_r(1), &
!                                     s_E(1),s_E(2),drytolerance,g,sw,fw)
            
            ! Eliminate ghost fluxes for wall
            do mw=1,3
                sw(mw)=sw(mw)*wall(mw)
                do m=1,3
                   fw(m,mw)=fw(m,mw)*wall(mw)
                enddo
            enddo

            ! Update speeds and waves
            ! Note that we represent all the waves in the first three arrays
            ! so it does not directly correspond to the two-layer case's wave
            ! structure
            s(i,:) = 0.d0
            fwave(i,:,:) = 0.d0
            
            s(i,1:3) = sw(:)
            fwave(i,1,1:3) = fw(1,:) * rho(1)
            fwave(i,n_index,1:3) = fw(2,:) * rho(1)
            fwave(i,t_index,1:3) = fw(3,:) * rho(1)
            
            ! Go on to next cell, lat-long and fluctuation calculations are 
            ! outside of this loop
            cycle
        ! ====================================================================
        ! Dry state, bottom layer to right
        else if(dry_state_r(2).and.(.not.dry_state_l(2))) then
            ! Assume currently no wetting or drying, otherwise we need to
            ! to check here if the state to the right is dry but shold be wet
            h_r(2) = h_l(2)
            hu_r(2) = -hu_l(2)
            u_r(2) = -u_l(2)
            hv_r(2) = hv_l(2)
            v_r(2) = v_l(2)
            
            flux_transfer_r = 0.d0
            flux_transfer_l = 0.d0
            momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r - h_l(2) - b_l)
            momentum_transfer(2) = 0.d0
            
        ! ====================================================================
        ! Dry state, bottom layer to left
        else if(dry_state_l(2).and.(.not.dry_state_r(2))) then
            h_l(2) = h_r(2)
            hu_l(2) = -hu_r(2)
            u_l(2) = -u_r(2)
            hv_l(2) = hv_r(2)
            v_l(2) = v_r(2)
            
            flux_transfer_r = 0.d0
            flux_transfer_l = 0.d0
            momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r + h_r(2) - b_l)
            momentum_transfer(2) = 0.d0
            
        ! ====================================================================
        ! Full two layer case
        else
            momentum_transfer(1) =  g * rho(1) * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
            momentum_transfer(2) = -g * rho(1) * h_ave(1) * (h_r(2) - h_l(2)) + g * rho(2) * h_ave(2) * (b_r - b_l)
            flux_transfer_r = g * rho(1) * h_r(1) * h_r(2)
            flux_transfer_l = g * rho(1) * h_l(1) * h_l(2)
        endif

        ! Check Richardson number
        kappa_l = (u_l(1) - u_l(2))**2 / (g*one_minus_r*sum(h_l))
        kappa_r = (u_r(1) - u_r(2))**2 / (g*one_minus_r*sum(h_r))
        if ((kappa_l > richardson_tolerance).and.(.not.dry_state_l(2))) then
            print "(a,i4,a,d16.8)","Hyperbolicity may have failed, kappa(",i,") = ",kappa_l
            print "(a,i4)","  Direction = ",ixy," Location = ",icom,jcom
        else if ((kappa_r > richardson_tolerance).and.(.not.dry_state_r(2))) then
            print "(a,i4,a,d16.8)","Hyperbolicity may have failed, kappa(",i,") = ",kappa_r
            print "(a,i4)","  Direction = ",ixy," Location = ",icom,jcom
        endif
        
        ! ====================================================================
        !  Calculate wave speeds
        ! ====================================================================
        ! Figure out which speed to use for transverse direction
        advected_speed(:) = 0.5d0 * (u_l(:) + u_r(:))
!         do j=1,2
!             if ((u_l(j) < 0.d0).and.(u_r(j) < 0.d0)) then
!                 advected_speed(j) = u_l(j)
!             else if ((u_l(j) > 0.d0).and.(u_r(j) > 0.d0)) then
!                 advected_speed(j) = u_r(j)
!             else
! !                 print *,"Transonic Riemann problem detected!"
! !                 print "(a,d16.8,a,d16.8,a)","(u_l,u_r) = (",u_l(j),",",u_r(j),")"
!                 advected_speed(j) = 0.5d0 * (u_l(j) + u_r(j))
!             endif
!         enddo
        if ((eigen_method > 0).and.(eigen_method <= 3)) then
            if (eigen_method == 1) then
                gamma_l = h_hat_l(2) / h_hat_l(1)
                gamma_r = h_hat_r(2) / h_hat_r(1)
                
                alpha(1) = 0.5d0*(gamma_l-1.d0+sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
                alpha(2) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
                alpha(3) = 0.5d0*(gamma_r-1.d0-sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
                alpha(4) = 0.5d0*(gamma_r-1.d0+sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
                
                s(i,1) = -sqrt(g*h_hat_l(1)*(1+alpha(1)))
                s(i,2) = -sqrt(g*h_hat_l(1)*(1+alpha(2)))
                s(i,5) = sqrt(g*h_hat_r(1)*(1+alpha(3)))
                s(i,6) = sqrt(g*h_hat_r(1)*(1+alpha(4)))
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
                stop "Eigenstructure calculation not implemented!"
            endif
            
            s(i,3:4) = advected_speed

            ! Compute eigenspace exactly based on eigenvalues provided
            eig_vec(1,:) = [1.d0,1.d0,0.d0,0.d0,1.d0,1.d0]
            
            eig_vec(n_index,:) = [s(i,1),s(i,2),0.d0,0.d0,s(i,5),s(i,6)]
            
            eig_vec(t_index,:) = [v_l(1),v_l(1),1.d0,0.d0,v_r(1),v_r(1)]

            eig_vec(4,1:2) = alpha(1:2)
            eig_vec(4,3:4) = 0.d0
            eig_vec(4,5:6) = alpha(3:4)
            
            eig_vec(n_index+3,:) = s(i,:) * eig_vec(4,:)
            
            eig_vec(t_index+3,1:2) = v_l(2) * alpha(1:2)
            eig_vec(t_index+3,3:4) = [0.d0,1.d0]
            eig_vec(t_index+3,5:6) = v_r(2) * alpha(3:4)
            
        else if (eigen_method == 4) then
            ! Directly compute eigenstructure using LAPACK if need be
            stop "Eigenstructure calculation not implemented!"            
        else if (eigen_method == 5) then
            ! Old eigenspeed and vector method with new alphas
            stop "Eigenstructure calculation not implemented!"
        else
            print "(a,i2,a)","Eigenstructure method ",eigen_method, &
                  " requested not available."
            stop 
        endif

        ! ====================================================================
        ! Compute jump in fluxes, this is generic due to terms set above
        do j=1,2
            layer_index = 3*(j-1)
            flux_r(layer_index+1) = rho(j) * hu_r(j)
            flux_r(layer_index+n_index) = rho(j) * (h_r(j) * u_r(j)**2 + 0.5d0 * g * h_r(j)**2)
            flux_r(layer_index+t_index) = rho(j) * h_r(j) * u_r(j) * v_r(j)
            
            flux_l(layer_index+1) = rho(j) * hu_l(j)
            flux_l(layer_index+n_index) = rho(j) * (h_l(j) * u_l(j)**2 + 0.5d0 * g * h_l(j)**2)
            flux_l(layer_index+t_index) = rho(j) * h_l(j) * u_l(j) * v_l(j)
        enddo
        ! Add extra flux terms
        flux_r(3 + n_index) = flux_r(3 + n_index) + flux_transfer_r
        flux_l(3 + n_index) = flux_l(3 + n_index) + flux_transfer_l
        
        delta = flux_r - flux_l
            
        ! Momentum transfer and bathy terms
        delta(n_index) = delta(n_index) + momentum_transfer(1)
        delta(n_index+3) = delta(n_index+3) + momentum_transfer(2)
        
        ! ====================================================================
        ! Project jump in fluxes - Use LAPACK's dgesv routine
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
            if (dry_state_l(2)) then
                print *,"left dry"
            else if (dry_state_r(2)) then
                print *,"right dry"
            endif
            print *,h_r(2),h_l(2)
            print *,hu_r(2),hu_l(2)
            print *,hv_r(2),hv_l(2)
            
            print "(a,i2)","In normal solver: ixy=",ixy
            print "(a,i3)","  Error solving R beta = delta,",info
            print "(a,i3,a,i3)","  Location: ",icom," ",jcom
            print "(a,6d16.8)","  Eigenspeeds: ",s(i,:)
            print "(a)","  Eigenvectors:"
            do j=1,6
                print "(a,6d16.8)","  ",(eig_vec(j,mw),mw=1,6)
            enddo
            stop
        endif
        beta = delta

        ! ====================================================================
        ! Compute fwaves
        forall(mw=1:mwaves)
            fwave(i,:,mw) = eig_vec(:,mw) * beta(mw)
        end forall
        
        ! ====================================================================
        !  Zero out waves that should be zero
        if (dry_state_r(2).and.(.not.dry_state_l(2))) then
!             do mw=1,6
!                 if (s(i,mw) > 0.d0) then
!                     do m=4,6
!                         if (fwave(i,m,mw) /= 0.d0) then
!                             print "(d16.8)",s(i,mw)
!                             print "(3d16.8)",(fwave(i,j,mw),j=4,6)
!                             fwave(i,4:6,mw) = 0.d0
!                         endif
!                     enddo
!                 endif
!             enddo
!             fwave(i,4:6,4) = 0.d0
!             print *,"======================================="
!             print *,"Right dry state...",i,jcom
!             do j=1,6
!                 print "(6d16.8)",(fwave(i,j,m),m=1,6)
!             enddo
!             fwave(i,4:6,4:6) = 0.d0
!             fwave(i,4:6,6) = 0.d0
!             print *,"======================================="
!             print *,"Right dry state..."
!             do j=1,6
!                 print "(6d16.8)",(fwave(i,j,m),m=1,6)
!             enddo
        else if (dry_state_l(2).and.(.not.dry_state_r(2))) then
!             fwave(i,4:6,1) = 0.d0
!             fwave(i,4:6,3) = 0.d0
!             print *,"======================================="
!             print *,"Left dry state..."
!             do j=1,6
!                 print "(6d16.8)",(fwave(i,j,m),m=1,6)
!             enddo
        endif
            
    enddo
    ! == End of Riemann Solver Loop per grid cell ============================
    
    ! ========================================================================
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
            h_r(1) = ql(i,1) / rho(1)
            h_l(1) = qr(i-1,1) / rho(1)
            h_r(2) = ql(i,4) / rho(2)
            h_l(2) = qr(i-1,4) / rho(2)
            if (h_r(2) < drytolerance .and. (.not.(h_l(2) < drytolerance))) then
                do m=4,6
                    if (apdq(i,m) /= 0.d0) then
                        print *,"========================"
                        print *,"Wave ",mw," equation ",m
                        print *,"s = ",s(i,mw)
                        print *,"f = ",fwave(i,m,mw)
                        print *,"amdq = ",(amdq(i,m))
                        print *,"apdq = ",(apdq(i,m))
                        stop
                    endif
                enddo
                apdq(i,4:6) = 0.d0
            endif
        enddo
    enddo

end subroutine rpn2