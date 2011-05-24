! ============================================================================
!  Program:     /Users/mandli/src/local/shallow_water/2d/multi_layer
!  File:        multilayer_module
!  Created:     2010-10-12
!  Author:      Kyle Mandli
! ============================================================================
!      Copyright (C) 2010-10-12 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module multilayer_module

    implicit none
    
    ! Physical parameters
    integer :: layers
    double precision, allocatable :: rho(:)
    double precision :: r,one_minus_r
    
    ! Algorithm parameters
    integer :: eigen_method,inundation_method
    double precision :: richardson_tolerance
    double precision, allocatable :: wave_tol(:)
    logical :: dry_limit
    
    ! Initial layer depths
    double precision, allocatable :: eta(:)
    double precision :: epsilon,sigma,init_location(2),angle
    integer :: init_type,wave_family
    
    ! Simple bathy states
    integer :: bathy_type = 2
    double precision :: bathy_location, bathy_left, bathy_right
    
    ! Complex bathy layout
    double precision :: x0,x1,x2,shelf_depth,basin_depth,beach_slope
    ! Calculated parameters
    double precision :: shelf_slope,eta_int
    
contains

    subroutine set_multilayer_params(data_file)

        implicit none
        character(len=*), optional, intent(in) :: data_file
        
        integer :: ios
        
        double precision :: A,B,step_height
        
        ! Open file
        if (present(data_file)) then
            open(unit=13,file=data_file,iostat=ios,status='old',action='read',access='sequential')            
            if ( ios /= 0 ) then
                print *,'Error opening "',data_file,'" for reading.'
                stop
            endif
        else
            open(unit=13,file='multilayer.data',iostat=ios,status='old',action='read',access='sequential')            
            if ( ios /= 0 ) then
                print *,'Error opening "multilayer.data" for reading.'
                stop
            endif
        endif
        
        dry_limit = .true.
    
        ! Physics parameters
        read(13,"(i3)") layers
        allocate(rho(layers))
        allocate(eta(layers))
        allocate(wave_tol(layers))
        read(13,*) rho
        if (layers > 1) then
            r = rho(1) / rho(2)
            one_minus_r = 1.d0 - r
        else
            r = -1.d0
            one_minus_r = 0.d0
        endif
        read(13,*)
        
        ! Algorithmic parameters
        read(13,"(i1)") eigen_method
        read(13,"(i1)") inundation_method
        read(13,"(d16.8)") richardson_tolerance
        read(13,*) wave_tol
        read(13,*) dry_limit
        read(13,*)
        
        ! Initial conditions
        read(13,*) eta
        read(13,*) init_type
        read(13,*) epsilon
        if (1 <= init_type .and. init_type <= 3) then  
            read(13,*) init_location
            read(13,*) wave_family
            if(2 <= init_type .and. init_type <= 3) then
                read(13,*) angle
                read(13,*) sigma
            endif
        endif
        
        ! Bathymetry
        read(13,*) bathy_type
        if (bathy_type == 0) then
            ! Probably should check here that mtopo is nonzero
            continue
        else if (bathy_type == 1) then
            read(13,*) bathy_location
            read(13,*) bathy_left
            read(13,*) bathy_right
        else if (bathy_type == 2 .or. bathy_type == 3) then
            read(13,*) x0
            read(13,*) x1
            read(13,*) x2
            read(13,*) basin_depth
            read(13,*) shelf_depth
            read(13,*) beach_slope
            read(13,*) step_height
            
            ! Calculated values
            if (layers > 1) then
                A = basin_depth - eta(2) + 0.5d0 * step_height
                B = shelf_depth - eta(2) - 0.5d0 * step_height
                eta_int = (A*x1 - B*x0) / (A-B)
                shelf_slope = A / (x0 - eta_int)
            else
                eta_int = 0.5d0 * (x0 + x1)
                shelf_slope = (basin_depth - shelf_depth) / (x0 - x1)
            endif
        endif
        
        close(13)

    end subroutine set_multilayer_params

end module multilayer_module
