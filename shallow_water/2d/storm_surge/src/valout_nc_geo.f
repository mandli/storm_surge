c
c -----------------------------------------------------
c     Routine to write netcdf files in the classic format
!        #jj-2011.02.04
!        # Each file written by the fortran code has 
!        # Dimensions:
!        #           timedimension : UNLIMITED
!        #           meqn          : The number of equations
!        #           dimx_<gridno> : X dimension for grid number <gridno>
!        #           dimy_<gridno> : Y dimension for grid number <gridno>
!        # Variables:
!        #           timedimension : Stores the time of the frame
!        #           ngrids        : Number of grids in this frame
!        #           naux          : Number of Auxilary Variables
!        #           ndim          : Number of Dimensions in the frame
!        #           grid_<gridno> : A grid of (dimx,dimy,meqn)
!        # Attributes:
!        # (grid_<no>) gridno      : The number of this grid <grid_no>
!        #           level         : The AMR level
!        #           dim_names     : a list of dimensions [dimx,dimy]
!        #           dim<x,y>.low  : The lowest dimension value 
!        #           dim<x,y>.d    : The distance between grid points 
c -----------------------------------------------------
c
      subroutine valout (lst, lend, time, nvar, naux)
c

      use multilayer_module, only: layers,rho
      
      implicit double precision (a-h,o-z)
      character*10  matname1, matname2

      include  "call.i"
      include 'netcdf.inc'
      real(kind=8) time
      integer ncid,rcode
      integer timeid,tVarID,meqnID,ngridsVarID,nauxVarID,ndimVarID
      integer dimxid,dimyid,xlowid,ylowid,dxid,dyid
      integer gridid
      integer ntimes
      character*2 gridstr
      character*40 dim_names
      REAL, ALLOCATABLE  ::grid(:,:,:)
      real dx,dy,xlow,ylow
      
      double precision :: h(2),hu(2),hv(2),eta(2),wind_x,wind_y
      
      iadd(i,j,ivar) = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1)
      iaddaux(i,j,iaux) = locaux + i - 1 + mitot*((iaux-1)*mjtot+j-1)

      ntimes=1
c ::::::::::::::::::::::::::: VALOUT ::::::::::::::::::::::::::::::::::;
c graphics output of soln values for contour or surface plots.
c modified for GeoClaw to output the surface level along with q.
c    surface = q(i,j,1) + aux(i,j,1)
c :::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::::::::::;

c
c     ### MATLAB graphics output
c
      if (matlabout) then
c        ###  make the file names and open output files
         matname1 = 'fort.qxxxx'
         matname2 = 'fort.txxxx'
         matunit1 = 50
         matunit2 = 60
         nstp     = matlabu
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            matname1(ipos:ipos) = char(ichar('0') + idigit)
            matname2(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue

        !!!!Define netcdf file
         rcode=NF_CREATE(matname1//'.nc',NF_NOCLOBBER,ncid)
         if(rcode.ne.NF_NOERR) print *,'ERROR OPENING NETCDF FILE'
         rcode=NF_DEF_DIM(ncid,'timedimension',NF_UNLIMITED,timeid)
         rcode=NF_DEF_VAR(ncid,'timedimension',NF_FLOAT,1,timeid,tVarID)
C          rcode=NF_DEF_DIM(ncid,'meqn',nvar+1,meqnid)
C        This is custom for the multilayer code
         rcode=NF_DEF_DIM(ncid,'meqn',4*layers+2,meqnid)
         rcode=NF_DEF_VAR(ncid,'ngrids',NF_INT,0,0,ngridsVarID)
         rcode=NF_DEF_VAR(ncid,'naux',NF_INT,0,0,nauxVarID)
         rcode=NF_DEF_VAR(ncid,'ndim',NF_INT,0,0,ndimVarID)
         rcode=NF_ENDDEF(ncid)
!         rcode=NF_CLOSE(ncid)
         if(rcode.ne.NF_NOERR) print *,'ERROR  LEAVING DEFINE MODE'
         level = lst
         ngrids = 0
 65      if (level .gt. lfine) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              rcode=NF_REDEF(ncid)
              if(rcode.ne.NF_NOERR) print *,'ERROR  REDEFINE MODE'
              write(gridstr,67) mptr
              
67            format(I2.2)              
              rcode=NF_DEF_DIM(ncid,'dimx_'//trim(gridstr),nx,dimxid)
              if(rcode.ne.NF_NOERR) print *,'ERROR  DEFINE DIMS'
              rcode=NF_DEF_DIM(ncid,'dimy_'//trim(gridstr),ny,dimyid)
              if(rcode.ne.NF_NOERR) print *,'ERROR  DEFINE DIMS'
              

              rcode=NF_DEF_Var(ncid,'grid_'//trim(gridstr),NF_FLOAT,4,
     &              (/dimxid,dimyid,meqnid,timeid/),gridid)
               if(rcode.ne.NF_NOERR) print *,'ERROR  DEFINE VAR'
               
              rcode=NF_PUT_ATT_INT(ncid,gridid,'gridno',NF_INT,1,
     &              mptr)
     
              rcode=NF_PUT_ATT_INT(ncid,gridid,'level',NF_INT,1,level)
              
              dim_names="['dimx','dimy']"
              rcode=NF_PUT_ATT_TEXT(ncid,gridid,'dim_names',
     &         LEN_TRIM(dim_names),TRIM(dim_names))
     
              rcode=NF_PUT_ATT_REAL(ncid,gridid,'dimx.lower',NF_DOUBLE,
     &              1,xlow)     
              rcode=NF_PUT_ATT_REAL(ncid,gridid,'dimy.lower',NF_DOUBLE,
     &              1,ylow)
     
              dx=hxposs(level)
              dy=hyposs(level)
              rcode=NF_PUT_ATT_REAL(ncid,gridid,'dimx.d',NF_FLOAT,1,
     &          dx)
              rcode=NF_PUT_ATT_REAL(ncid,gridid,'dimy.d',NF_FLOAT,1,
     &          dy) 
     
              rcode=NF_ENDDEF(ncid)
              
C               allocate(grid(nx,ny,nvar+1))
C             Changed number of output in q grid
              allocate(grid(nx,ny,4*layers+2))
              
!!! with netcdf4 we can have groups that will encapsulate all of this data, but for now....use the netcdf3 format "universal" with scientific python

      do j = nghost+1, mjtot-nghost
         do i = nghost+1, mitot-nghost
            do ivar=1,nvar
               if (dabs(alloc(iadd(i,j,ivar))) .lt. 1d-90) then
                  alloc(iadd(i,j,ivar)) = 0.d0
               endif
               if (layers > 1) then
                   do k=1,layers
                       index = 3*(k-1)
                       h(k) = alloc(iadd(i,j,index+1)) / rho(k)
                       hu(k) = alloc(iadd(i,j,index+2)) / rho(k)
                       hv(k) = alloc(iadd(i,j,index+3)) / rho(k)
                   enddo
                   eta(2) = h(2) + alloc(iaddaux(i,j,1))
               else
                   k =1 
                   index = 3*(k-1)
                   h(k) = alloc(iadd(i,j,index+1))
                   hu(k) = alloc(iadd(i,j,index+2))
                   hv(k) = alloc(iadd(i,j,index+3))
                   eta(2) = alloc(iaddaux(i,j,1))
               endif
               eta(1) = h(1) + eta(2)
               wind_x = alloc(iaddaux(i,j,4))
               if (abs(wind_x) < 1d-90) then
                   wind_x = 0.d0
               endif
               wind_y = alloc(iaddaux(i,j,5))
               if (abs(wind_y) < 1d-90) then
                   wind_y = 0.d0
               endif
               
               do k=1,layers
                   index = 3*(k-1)
                   grid(i-nghost,j-nghost,index+1) = h(k)
                   grid(i-nghost,j-nghost,index+2) = hu(k)
                   grid(i-nghost,j-nghost,index+2) = hv(k)
               enddo
               do k=1,layers
                   index = 3*(k-1)
                   grid(i-nghost,j-nghost,index+3*layers) = eta(k)
               enddo
               grid(i-nghost,j-nghost,4*layers+1) = wind_x
               grid(i-nghost,j-nghost,4*layers+2) = wind_yl
            enddo
         enddo
      enddo
      
      rcode=NF_PUT_VARA_REAL(ncid,gridid,(/1,1,1,1/), 
     &      (/nx,ny,nvar+1,1/),grid)

      deallocate(grid)

109       format(4e26.16)
            mptr = node(levelptr, mptr)
            go to 70
 80      level = level + 1
         go to 65

 90     continue

      rcode=NF_PUT_VAR_DOUBLE(ncid,tVarID,time)
      if(rcode.ne.NF_NOERR) print *,'ERROR  Write Time'
      rcode=NF_PUT_VAR_INT(ncid,ngridsVarID,int(ngrids))
      if(rcode.ne.NF_NOERR) print *,'ERROR  Write GridNo'
      rcode=NF_PUT_VAR_INT(ncid,nauxVarID,3)
      rcode=NF_PUT_VAR_INT(ncid,ndimVarID,2)
      rcode=NF_CLOSE(ncid)
      

      write(6,601) matlabu,time
  601 format('GeoClaw: Frame ',i4,
     &       ' output files done at time t = ', d12.6,/)

      matlabu = matlabu + 1

      endif

      return
      end

