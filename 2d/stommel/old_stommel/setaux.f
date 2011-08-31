c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # Set the aux array:
c
c     # aux(i,j,1) = average normal velocity along left edge
c     # aux(i,j,2) = average normal velocity along bottom edge
c     # aux(i,j,3) = capacity of cell: fraction of cell in flow domain
c     #              with some fixes for small triangular cells
c
      implicit real*8(a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 3)
c 
c     # local storage  assumes mbc = 2, maxmx and maxmy < maxm2 - 2
      parameter (maxmbc = 6)
      parameter (maxm2 = 200)
      dimension frleft(1-maxmbc:maxm2+maxmbc,1-maxmbc:maxm2+maxmbc)
      dimension frbot(1-maxmbc:maxm2+maxmbc,1-maxmbc:maxm2+maxmbc)
c
c
      if (mx .gt. maxm2 .or. my .gt. maxm2) then
	 write(6,*) '*** ERROR in setaux: mx,my,maxm2:',mx,my,maxm2
	 stop
	 endif
      if (mbc .gt. maxmbc) then
	 write(6,*) '*** ERROR in setaux: mbc,maxmbc:',mbc,maxmbc
	 stop
	 endif

      do 20 j=1-mbc,my+mbc
         do 20 i=1-mbc,mx+mbc
c
c           # coordinates of lower-left corner and adjacent corners:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy
            xup = xll + dx
            yup = yll + dy
c
c           # compute average normal velocity across each edge by
c           # differencing stream function:
c
            aux(i,j,1) = (stream(xll,yll) - stream(xll,yup)) / dy
            aux(i,j,2) = (stream(xup,yll) - stream(xll,yll)) / dx
c
c           # compute fraction of cell area, left edge, bottom edge,
c           # interior to physical domain
c
            call cartfr(xll,yll,dx,dy,frarea,frleft(i,j),frbot(i,j))
c           if (frarea.lt.0.1d0) frarea = 0.1d0
            aux(i,j,3) = frarea
c	    write(17,1701) i,j,aux(i,j,1),aux(i,j,2),aux(i,j,3)
 1701       format(2i4,3d16.6)
c
   20       continue

c
c     # triangular cell fix-up:
c
      do 30 j=1-mbc,my+mbc-1
         do 30 i=1-mbc,mx+mbc-1
            hmax = dmax1(frleft(i,j),frbot(i,j),
     &                   frleft(i+1,j),frbot(i,j+1))
            if (hmax .gt. 0.d0) then
                aux(i,j,3) = aux(i,j,3) / hmax
              endif
   30       continue
c
      do  j=1-mbc,my+mbc
         do  i=1-mbc,mx+mbc
c	    aux(i,j,3) = dmax1(aux(i,j,3), 1d-3)
            aux(i,j,3) = 1.d0
c	    write(12,*) i,j,aux(i,j,3)
	    enddo
	    enddo
c
      return
      end

