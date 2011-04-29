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
    integer :: eigen_method
    double precision :: richardson_tolerance
    double precision, allocatable :: wave_tol(:)
    logical :: dry_limit
    
    ! Initial layer depths
    double precision, allocatable :: eta(:)
    double precision :: epsilon,sigma,init_location(2)
    integer :: init_type,wave_family
    
    ! Simple bathy states
    double precision :: bathy_location, bathy_left, bathy_right
    
contains

    subroutine set_multilayer_params(data_file)

        implicit none
        character(len=*), optional, intent(in) :: data_file
        
        integer :: ios
        
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
        read(13,*) layers
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
        read(13,"(i2)") eigen_method
        read(13,"(d16.8)") richardson_tolerance
        read(13,*) wave_tol
        read(13,*) dry_limit
        read(13,*)
        
        ! Initial conditions
        read(13,*) eta
        read(13,*) init_type
        read(13,*) init_location
        read(13,*) wave_family
        read(13,*) epsilon
        read(13,*) sigma
        
        ! Bathymetry
        read(13,*) bathy_location
        read(13,*) bathy_left
        read(13,*) bathy_right
        
        close(13)

    end subroutine set_multilayer_params

end module multilayer_module
