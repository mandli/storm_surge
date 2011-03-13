subroutine b4step1(maxmx,mbc,mx,meqn,q,xlower,dx,t,dt,maux,aux)

    use parameters_module

    implicit none
    
    ! Arguments
    integer, intent(in) :: maxmx,mbc,mx,meqn,maux
    double precision, intent(in) :: xlower,dx,t,dt
    double precision, intent(inout) :: q(1-mbc:maxmx+mbc,meqn)
    double precision, intent(inout) :: aux(1-mbc:maxmx+mbc,maux)
    
    ! Calculate wind field
    call set_wind(maxmx,mbc,mx,xlower,dx,t,aux(:,2))
    
end subroutine b4step1