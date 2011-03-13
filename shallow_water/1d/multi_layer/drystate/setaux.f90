! ============================================================================
!  Sets aux arrays
!    aux(:,1) = Bathymetry
!    aux(:,2) = Initial steady state, top layer depth, h_1 = eta_1 - eta_2
!    aux(:,3) = Initial steady state, bottom layer depth, h_2 = eta_2 - b
! ============================================================================
subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
    
    use parameters_module
    
    implicit none
    
    ! Input parameters
    integer, intent(in) :: maxmx,mbc,mx,maux
    double precision, intent(in) :: xlower,dx
    
    ! Auxillary array  
    double precision, intent(inout) :: aux(1-mbc:maxmx+mbc, maux)
     
    ! Local
    integer, parameter :: out_unit = 13
    integer :: i, ios
    double precision :: x
    
    do i=1-mbc,mx+2*mbc
        x = xlower+(i-0.5)*dx
        if (x < bathy_location) then
            aux(i,1) = bathy_left
        else
            aux(i,1) = bathy_right
        endif
    enddo
    
    ! Wind
    aux(:,2) = 0.d0
    
    ! Write out auxillary array
    print *,"Outputting bathymetry to 'fort.aux' file."
    open(unit=out_unit, file='fort.aux', iostat=ios, action="write")
    if ( ios /= 0 ) stop "Error opening file name"
    
    do i=1,mx
        write(unit=out_unit, fmt="(d16.8)") aux(i,1)
    enddo
    
    close(out_unit)
    
end subroutine setaux
