! ============================================================================
!  Sets aux arrays
!    aux(:,1) = Bathymetry
!    aux(:,2) = Wind field
!    aux(:,3) = Initial steady state, top layer depth, h_1 = eta_1 - eta_2
!    aux(:,4) = Initial steady state, bottom layer depth, h_2 = eta_2 - b
! ============================================================================
subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
    
    use parameters_module
    
    implicit none
    
    ! Input parameters
    integer, intent(in) :: maxmx,mbc,mx,maux
    double precision, intent(in) :: xlower,dx
    
    ! Auxillary array  
    double precision, intent(inout) :: aux(maux,1-mbc:maxmx+mbc)
     
    ! Local
    integer, parameter :: out_unit = 13
    integer :: i, ios
    double precision :: x
    
    ! maux hack
    ml_maux = maux
    
    do i=1-mbc,mx+mbc
        x = xlower+(i-0.5)*dx
        ! Jump in bathymetry
        if (bathy_type == 1) then
            if (x < bathy_location) then
                aux(1,i) = bathy_left
            else
                aux(1,i) = bathy_right
            endif
        ! Simple shelf
        else if (bathy_type == 2) then
            if (x < x0) then
                aux(1,i) = basin_depth
            else if (x0 <= x .and. x < x1) then
                aux(1,i) = shelf_slope * (x-x0) + basin_depth
            else if (x1 <= x) then
                aux(1,i) = shelf_depth
            endif
        endif
        
        ! Set initial states
        if (x < init_location) then
            if (eta_left(2) > aux(1,i)) then
                aux(3,i) = eta_left(1) - eta_left(2)
                aux(4,i) = eta_left(2) - aux(1,i)
            else
                aux(3,i) = eta_left(1) - aux(1,i)
            endif
        else
            if (eta_right(2) > aux(1,i)) then
                aux(3,i) = eta_right(1) - eta_right(2)
                aux(4,i) = eta_right(2) - aux(1,i)
            else
                aux(3,i) = eta_right(1) - aux(1,i)
            endif
        endif
    enddo
    
    
    
    ! Calculate initial wind field
    call set_wind(maxmx,mbc,mx,xlower,dx,0.d0,aux)
    
    ! Write out auxillary array
    print *,"Outputting bathymetry to 'fort.aux' file."
    open(unit=out_unit, file='fort.aux', iostat=ios, action="write")
    if ( ios /= 0 ) stop "Error opening file name"
    
    do i=1,mx
        write(unit=out_unit, fmt="(d16.8)") aux(1,i)
    enddo
    
    close(out_unit)
    
end subroutine setaux
