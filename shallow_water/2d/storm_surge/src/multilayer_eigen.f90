! ============================================================================
!  File:        multilayer_eigen.f90
!  Created:     2011-05-04
!  Author:      Kyle Mandli
! ============================================================================
!      Copyright (C) 2011-05-04 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

subroutine linearized_eigen(h_l,h_r,hu_l,hu_r,hv_l,hv_r,u_l,u_r,v_l,v_r, &
                            n_index,t_index,s,eig_vec)

    use multilayer_module, only: r
    use geoclaw_module, only: grav

    implicit none
    
    ! Input
    double precision, dimension(2), intent(in) :: h_l,h_r,hu_l,hu_r,hv_l,hv_r
    double precision, dimension(2), intent(in) :: u_l,u_r,v_l,v_r
    integer, intent(in) :: n_index,t_index
    
    ! Output
    double precision, intent(inout) :: s(6),eig_vec(6,6)
        
    ! Local
    double precision :: gamma_l,gamma_r,alpha(4),g
    
    g = grav
        
    ! Calculate relevant quantities
    gamma_l = h_l(2) / h_l(1)
    gamma_r = h_r(2) / h_r(1)

    alpha(1) = 0.5d0*(gamma_l-1.d0+sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
    alpha(2) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
    alpha(3) = 0.5d0*(gamma_r-1.d0-sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
    alpha(4) = 0.5d0*(gamma_r-1.d0+sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))

    s(1) = -sqrt(g*h_l(1)*(1+alpha(1)))
    s(2) = -sqrt(g*h_l(1)*(1+alpha(2)))
    s(3:4) = 0.5d0 * (u_l(:) + u_r(:))
    s(5) = sqrt(g*h_r(1)*(1+alpha(3)))
    s(6) = sqrt(g*h_r(1)*(1+alpha(4)))

    ! Compute eigenspace exactly based on eigenvalues provided
    eig_vec(1,:) = [1.d0,1.d0,0.d0,0.d0,1.d0,1.d0]
    
    eig_vec(n_index,:) = [s(1),s(2),0.d0,0.d0,s(5),s(6)]
    
    eig_vec(t_index,:) = [v_l(1),v_l(1),1.d0,0.d0,v_r(1),v_r(1)]

    eig_vec(4,1:2) = alpha(1:2)
    eig_vec(4,3:4) = 0.d0
    eig_vec(4,5:6) = alpha(3:4)
    
    eig_vec(n_index+3,:) = s * eig_vec(4,:)
    
    eig_vec(t_index+3,1:2) = v_l(2) * alpha(1:2)
    eig_vec(t_index+3,3:4) = [0.d0,1.d0]
    eig_vec(t_index+3,5:6) = v_r(2) * alpha(3:4)

end subroutine linearized_eigen

subroutine vel_diff_eigen(h_l,h_r,hu_l,hu_r,hv_l,hv_r,u_l,u_r,v_l,v_r, &
                            n_index,t_index,s,eig_vec)

    use multilayer_module, only: one_minus_r,r
    use geoclaw_module, only: grav

    implicit none
    
    ! Input
    double precision, dimension(2), intent(in) :: h_l,h_r,hu_l,hu_r,hv_l,hv_r
    double precision, dimension(2), intent(in) :: u_l,u_r,v_l,v_r
    integer, intent(in) :: n_index,t_index
    
    ! Output
    double precision, intent(inout) :: s(6),eig_vec(6,6)
        
    ! Local
    double precision :: total_depth_l,total_depth_r,mult_depth_l,mult_depth_r
    
    total_depth_l = sum(h_l)
    total_depth_r = sum(h_r)
    mult_depth_l = product(h_l)
    mult_depth_r = product(h_r)
                      
    s(1) = - sqrt(grav*total_depth_l) + 0.5d0 * mult_depth_l/ total_depth_l**(3/2) * one_minus_r
    s(2) = - sqrt(grav * mult_depth_l / total_depth_l * one_minus_r)
    s(3:4) = 0.5d0 * (u_l + u_r)
    s(5) = sqrt(grav * mult_depth_r / total_depth_r * one_minus_r)
    s(6) = sqrt(grav*total_depth_r) - 0.5d0 * mult_depth_r / total_depth_r**(3/2) * one_minus_r

    eig_vec(1,:) = [1.d0,1.d0,0.d0,0.d0,1.d0,1.d0]
    eig_vec(n_index,:) = [s(1),s(2),0.d0,0.d0,s(5),s(6)]
    eig_vec(t_index,:) = [v_l(1),v_l(1),1.d0,0.d0,v_r(1),v_r(1)]
    eig_vec(4,1:2) = ((s(1:2)-u_l(1))**2 - grav*h_l(1)) / (r*grav*h_l(1))
    eig_vec(4,3:4) = 0.d0
    eig_vec(4,5:6) = ((s(5:6)-u_r(1))**2 - grav*h_r(1)) / (r*grav*h_r(1))
    eig_vec(n_index+3,:) = s(:) * eig_vec(4,:)
    eig_vec(t_index+3,1:2) = v_l(2) * eig_vec(4,1:2)
    eig_vec(t_index+3,3:4) = [0.d0,1.d0]
    eig_vec(t_index+3,5:6) = v_r(2) * eig_vec(4,5:6)
    
end subroutine vel_diff_eigen


subroutine lapack_eigen(h_l,h_r,hu_l,hu_r,hv_l,hv_r,u_l,u_r,v_l,v_r, &
                            n_index,t_index,s,eig_vec)

    use multilayer_module, only: r
    use geoclaw_module, only: drytolerance,grav

    implicit none
    
    ! Input
    double precision, dimension(2), intent(in) :: h_l,h_r,hu_l,hu_r,hv_l,hv_r
    double precision, dimension(2), intent(in) :: u_l,u_r,v_l,v_r
    integer, intent(in) :: n_index,t_index
    
    ! Output
    double precision, intent(inout) :: s(6),eig_vec(6,6)
    
    ! Local
    integer, parameter :: lwork = 6*6
    integer :: i,j,m,info
    double precision :: h_ave(2),u_ave(2),v_ave(2),A(6,6),A_copy(6,6)
    double precision :: imag_evalues(6),empty,work(1,lwork)
    double precision :: g
    
    g = grav

    ! Construct flux matrix
    A = 0.d0
    h_ave(:) = 0.5d0 * (h_l(:) + h_r(:))
    u_ave(:) = 0.5d0 * (u_l(:) + u_r(:))
    v_ave(:) = 0.5d0 * (v_l(:) + v_r(:))
    ! We need to do this since one of the rows cannot be swapped
    if (n_index == 2) then
        A(t_index,1:3) = [-u_ave(1)*v_ave(1),v_ave(1),u_ave(1)]
        A(t_index+3,4:6) = [-u_ave(2)*v_ave(2),v_ave(2),u_ave(2)]
    else
        A(t_index,1:3) = [-u_ave(1)*v_ave(1),u_ave(1),v_ave(1)]
        A(t_index+3,4:6) = [-u_ave(2)*v_ave(2),u_ave(2),v_ave(2)]
    endif
        
    A(1,n_index) = 1.d0
    
    A(n_index,1) = -u_ave(1)**2 + g*h_ave(1)
    A(n_index,n_index) = 2*u_ave(1)
    A(n_index,4) = r*g*h_ave(1)
    
    A(4,n_index+3) = 1.d0
    
    A(n_index+3,1) = g*h_ave(2)
    A(n_index+3,4) = -u_ave(2)**2 + g*h_ave(2)
    A(n_index+3,n_index+3) = 2.d0 * u_ave(2)
    
    A_copy = A
    
    ! Call LAPACK
    call dgeev('N','V',6,A,6,s,imag_evalues,empty,1,eig_vec,6,work,lwork,info)
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
    do i=1,6
        if (eig_vec(1,i) /= 0.d0) then
            eig_vec(:,i) = eig_vec(:,i) / eig_vec(1,i)
        endif
        if (abs(imag_evalues(i)) > 0.d0) then
            print "(a,i1,a,d16.8)","Imaginary eigenvalue(",i,") > 0.0",imag_evalues(i)
            stop
        endif
    enddo

end subroutine lapack_eigen

! ============================================================================
!  Single layer eigensolver
!   Note that this routine puts the result in the upper 3x3 matrix, the rest
!   is left as zeros.
! ============================================================================
subroutine single_layer_eigen(h_l,h_r,hu_l,hu_r,hv_l,hv_r,u_l,u_r,v_l,v_r, &
                            n_index,t_index,s,eig_vec)

    use geoclaw_module, only: grav

    implicit none
    
    ! Input
    double precision, dimension(2), intent(in) :: h_l,h_r,hu_l,hu_r,hv_l,hv_r
    double precision, dimension(2), intent(in) :: u_l,u_r,v_l,v_r
    integer, intent(in) :: n_index,t_index
    
    ! Output
    double precision, intent(inout) :: s(6),eig_vec(6,6)

    s = 0.d0
    s(1) = u_l(1) - sqrt(grav*h_l(1))
    s(2) = 0.5d0 * (u_r(1) + u_l(1))
    s(3) = u_r(1) + sqrt(grav*h_r(1))
    
    eig_vec = 0.d0
    eig_vec(1,1:3) = [1.d0,0.d0,1.d0]
    eig_vec(n_index,1:3) = [s(1),0.d0,s(3)]
    eig_vec(t_index,1:3) = [v_l(1),1.d0,v_r(1)]
    
end subroutine single_layer_eigen