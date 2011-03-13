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
     
    ! Constant bathy
    select case(abs(bathy_type))
        ! Constant
        case(1)
            aux(:,1) = -sign(1,bathy_type) * depth
        case(2)
            forall(i=1:mx,xlower+(i-0.5d0)*dx < init_location)
                aux(i,1) = -1.0d0
            end forall
            forall(i=1:mx,xlower+(i-0.5d0)*dx >= init_location)
                aux(i,1) = -0.3d0
            end forall
            
        case default
            print *,"Invalid bathymetry type ", bathy_type
    end select
            

    ! Calculate initial wind field
    call set_wind(maxmx,mbc,mx,xlower,dx,-99.d99,aux(:,2))
    
    ! Write out auxillary array
    print *,"Outputting bathymetry to 'fort.aux' file."
    open(unit=out_unit, file='fort.aux', iostat=ios, action="write")
    if ( ios /= 0 ) stop "Error opening file name"
    
    do i=1,mx
        write(unit=out_unit, fmt="(d16.8)") aux(i,1)
    enddo
    
    close(out_unit)
    
end subroutine setaux
