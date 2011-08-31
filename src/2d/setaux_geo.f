c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlow,ylow,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays
c
c     # aux(i,j,1) = Z(x,y) topography
c     #                     (negative below sea level for topoymetry)
c
c     # If icoordsys=2 then lat-lon coordinates on the sphere and
c     #    aux(i,j,2) = area ratio (capacity function -- set mcapa = 2)
c     #    aux(i,j,3) = length ratio for edge
c
c     # Storm surge related auxillary variables
c     # aux(i,j,4:5) = wind speed in the x and y directions
c     # aux(i,j,6) = pressure
c     # aux(i,j,7:8) = Initial depths


      use multilayer_module
      use hurricane_module
      use geoclaw_module
      use topo_module
      
      implicit double precision (a-h,o-z)

      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      include "call.i"

      if (icoordsys.eq.2) then
         if (mcapa .ne. 2 .or. maux.lt.3) then
            write(6,*) 'ERROR in setaux:  for icoordsys=2'
            write(6,*) '     need mcapa = 2 and maux >= 3'
            write(6,*) '     have mcapa = ',mcapa,'  maux = ',maux
            stop
            endif
         endif

      do j=1-mbc,my+mbc
         ycell = ylow +(j-0.5d0)*dy
         yjm = ylow +(j-1.d0)*dy
         yjp = ylow + j*dy

         do i=1-mbc,mx+mbc
            xcell= xlow + (i- 0.5d0)*dx
            xim = xlow + (i - 1.d0)*dx
            xip = xlow + i*dx

            if (icoordsys.eq.2) then
c           # for lat-lon grid on sphere:
               deg2rad = pi/180.d0
               aux(i,j,2)= deg2rad*Rearth**2*
     &               (sin(yjp*deg2rad)-sin(yjm*deg2rad))/dy
               aux(i,j,3)= yjm*deg2rad
            else
               aux(i,j,2) = 1.d0
               aux(i,j,3) = 1.d0
	            endif
	        
            if (mtopofiles.gt.0) then
               topoint=0.d0
               call cellgridintegrate(topoint,xim,xcell,xip,yjm,ycell,
     &	        yjp,xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo,
     &	        mxtopo,mytopo,mtopo,i0topo,mtopoorder,
     &	        mtopofiles,mtoposize,topowork)
               aux(i,j,1) = topoint/(dx*dy*aux(i,j,2))
            else
                stop "No topography files found!"
            endif
            
            enddo
         enddo
         
         

c     Initialize wind and pressure auxillary variables
      call hurricane_wind(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                      -ramp_up_time,aux(:,:,4:5))
      call hurricane_pressure(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,
     &                         dy,-ramp_up_time,aux(:,:,6))

c     This actually only handles 2 layers right now...
      if (layers > 1) then
          do i=1-mbc,mx+mbc
              do j=1-mbc,my+mbc
                  if (eta(2) > aux(i,j,1)) then
                      aux(i,j,7) = eta(1) - eta(2)
                      aux(i,j,8) = eta(2) - aux(i,j,1)
                  else
                      aux(i,j,7) = eta(1) - aux(i,j,1)
                      aux(i,j,8) = 0.d0
                  endif
              enddo
          enddo
      endif
      
      return

c     -----------------------------------------------------------------
c     # output aux array for debugging:
      open(23, file='fort.aux',status='unknown',form='formatted')
C       write(23,*) 'Setting aux arrays'
C       write(23,*) ' '
      j = 10
      do i=1,mx
           print *,i,j,(aux(i,j,m),m=1,1)
           write(23,*) i,j,(aux(i,j,m),m=1,1)
        enddo
 231  format(2i4,4d16.8)
      close(23)
      stop
      return
      end