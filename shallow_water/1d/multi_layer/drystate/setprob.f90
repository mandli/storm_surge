module parameters_module

    implicit none
    
    ! Physics constants
    double precision, parameter :: PI = 3.141592654d0, g = 9.8d0
    double precision :: rho(2)
    double precision :: r, one_minus_r
    
    ! Algorithm settings
    double precision :: dry_tolerance

    ! Initial condition settings
    integer :: wave_family
    double precision :: init_location,eta_1,eta_2,epsilon
    double precision :: bathy_location,bathy_left,bathy_right

contains

    ! ========================================================================
    !  Setup all parameters in the module
    ! ========================================================================
    subroutine setup_parameters(data_file)

        implicit none

        ! Input
        character(len=*), optional, intent(in) :: data_file

        ! Locals
        integer :: ios,i

        ! Open file
        if (present(data_file)) then
            call opendatafile(13,data_file)
        else
            call opendatafile(13,'problem.data')
        endif

        ! Physics constants
        read(13,*) rho(1)
        read(13,*) rho(2)
        r = rho(1) / rho(2)
        one_minus_r = 1.d0 - r
        read(13,*)
        
        ! Algorithm settings
        read(13,*) dry_tolerance
        read(13,*)
        
        ! Initial conditions
        read(13,*) init_location
        read(13,*) wave_family
        read(13,*) eta_1
        read(13,*) eta_2
        read(13,*) epsilon
        read(13,*)
        
        ! Bathymetry
        read(13,*) bathy_location
        read(13,*) bathy_left
        read(13,*) bathy_right
        
    end subroutine setup_parameters

end module parameters_module


subroutine setprob()

    use parameters_module

    implicit none
    
    call setup_parameters()

end subroutine setprob