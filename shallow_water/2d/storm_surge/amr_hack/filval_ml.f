c
c ------------------------------------------------------------------
c
      subroutine filval(val,mitot,mjtot,hx,hy,lev,time,
     1                  valc,auxc,mic,mjc,
     2                  xleft,xright,ybot,ytop,nvar,
     3                  mptr,ilo,ihi,jlo,jhi,aux,naux,locflip)
 
      use geoclaw_module
 
      implicit double precision (a-h,o-z)

      include "call.i"

      dimension   val(mitot,mjtot,nvar), valc(mic,mjc,nvar)
      dimension   aux(mitot,mjtot,naux), auxc(mic,mjc,naux)
      dimension   dudx(max1d), dudy(max1d)

      dimension coarseq(3,3,3)
      dimension coarseeta(3,3)
      dimension coarseb(3,3)
      dimension coarseinterp(3,3,3)
      dimension velfine(2)
      dimension coarsevelmax(2)
      dimension coarsevelmin(2)

      logical redo

c
c :::::::::::::::::::::::::::::: FILVAL ::::::::::::::::::::::::::
c
c create and fill coarser (lev-1) patch with one extra coarse cell all
c around, plus the ghost cells . will interpolate from this patch to grid mptr 
c without needing special boundary code. 
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     # indext into eta array for surface values:
      iaddeta(i,j) = loceta + i-1 + mic*(j-1)

      levc    = lev - 1
      lratiox = intratx(levc)
      lratioy = intraty(levc)
      hxcrse  = hx*lratiox
      hycrse  = hy*lratioy
      xl      = xleft  - hxcrse
      xr      = xright + hxcrse
      yb      = ybot   - hycrse
      yt      = ytop   + hycrse
c
c     set integer indices for coarser patch enlarged by 1 cell 
c     (can stick out of domain). proper nesting will insure this one
c     call is sufficient.
      iclo   = ilo/lratiox - 1
      jclo   = jlo/lratioy - 1
      ichi   = (ihi+1)/lratiox - 1 + 1
      jchi   = (jhi+1)/lratioy - 1 + 1
      ng     = 0

c    :::  mcapa  is the capacity function index

      if (naux .eq. 0) then
c     if (mcapa .eq. 0) then
        if (xperdom .or. yperdom .or. spheredom) then
          call preintcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,
     &                    locflip)
        else
          call intcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,1,1)
        endif
      else  ! intersect grids and copy all (soln and aux)
        if (xperdom .or. yperdom .or. spheredom) then
          call preicall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,
     &                  jchi,levc,locflip)
        else
          call icall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi,
     &               levc,1,1)
        endif
      endif
      call bc2amr(valc,auxc,mic,mjc,nvar,naux,
     1            hxcrse,hycrse,levc,time,
     2            xl,xr,yb,yt,
     3            xlower,ylower,xupper,yupper,
     4            xperdom,yperdom,spheredom)



c-----------------------------
c     # for shallow water over topography,
c     # in coarse cells convert from h,
c     # to eta,  before interpolating:
c-----------------------------
      toldry = drytolerance

c     #prepare slopes - use min-mod limiters
      do j=2, mjc-1
      do i=2, mic-1

       	redo=.false.

c     # coarse values surrounding the current coarse value are set
c     # 5 coarse values surrounding and including the current coarse value are set
         coarsevelmax(1)=-1.d99
         coarsevelmax(2)=-1.d99
         coarsevelmin(1)=1.d99
         coarsevelmin(2)=1.d99

         wetcount=0.d0
         etaave=0.d0
         sxlimit=1.d0
         sylimit=1.d0

         do ii=-1,1
         do jj=-1,1
         if (ii*jj.eq.0) then 
             coarseb(2+ii,2+jj)= auxc(i+ii,j+jj,1)
             do ivar=1,nvar
                   coarseq(2+ii,2+jj,ivar)= valc(i+ii,j+jj,ivar)                
             enddo
             coarseq(2+ii,2+jj,1)=max(0.d0,coarseq(2+ii,2+jj,1))
             coarseeta(2+ii,2+jj)=coarseq(2+ii,2+jj,1)
     &                          +coarseb(2+ii,2+jj)
             if (coarseq(2+ii,2+jj,1).le.toldry) then
                    coarseu=0.d0
                    coarsev=0.d0
                    coarseq(2+ii,2+jj,2)=0.d0
                    coarseq(2+ii,2+jj,3)=0.d0
                    coarseinterp(2+ii,2+jj,1)=sealevel
             else
                    wetcount=wetcount+1.d0
                    etaave=etaave + coarseeta(2+ii,2+jj)
                    coarseu=coarseq(2+ii,2+jj,2)/coarseq(2+ii,2+jj,1)
                    coarsev=coarseq(2+ii,2+jj,3)/coarseq(2+ii,2+jj,1)
                    coarseinterp(2+ii,2+jj,1)=coarseeta(2+ii,2+jj)
             endif

             coarseinterp(2+ii,2+jj,2)=coarseq(2+ii,2+jj,2)
             coarseinterp(2+ii,2+jj,3)=coarseq(2+ii,2+jj,3)

             coarsevelmax(1)=max(coarseu,coarsevelmax(1))
             coarsevelmin(1)=min(coarseu,coarsevelmin(1))
             coarsevelmax(2)=max(coarsev,coarsevelmax(2))
             coarsevelmin(2)=min(coarsev,coarsevelmin(2))
         endif
         enddo
         enddo


c:::::::::::::::::::limited slope interpolation ::::::::::::::::::::
         do ivar = 1,nvar
            s1p=coarseinterp(3,2,ivar)-coarseinterp(2,2,ivar)
            s1m=coarseinterp(2,2,ivar)-coarseinterp(1,2,ivar)
            slopex=dmin1(dabs(s1p),dabs(s1m))*dsign(1.d0,
     &          coarseinterp(3,2,ivar)-coarseinterp(1,2,ivar))
            if (s1m*s1p.le.0.d0) then
                slopex=0.d0
            endif

            s1p=coarseinterp(2,3,ivar)-coarseinterp(2,2,ivar)
            s1m=coarseinterp(2,2,ivar)-coarseinterp(2,1,ivar)
            slopey=dmin1(dabs(s1p),dabs(s1m))*dsign(1.d0,
     &                 coarseinterp(2,3,ivar)-coarseinterp(2,1,ivar))
            if (s1m*s1p.le.0.d0) then
                 slopey=0.d0
            endif

            if (.false.) then 
            do kk=-1,1
               if (coarseq(2+kk,2,1).lt.toldry) then
                  sxlimit=0.d0
               endif

               if (coarseq(2,2+kk,1).lt.toldry) then
                  sylimit=0.d0
               endif
             enddo
             endif

             dudx(i) = slopex
             dudy(i) = slopey


c     interp. from coarse cells to fine grid

             do ico = 1,lratiox
             do jco = 1,lratioy

                yoff  = (float(jco) - .5)/lratioy - .5
	          xoff = (float(ico) - .5)/lratiox - .5
                jfine = (j-2)*lratioy + nghost + jco
                ifine   = (i-2)*lratiox + nghost + ico
                val(ifine,jfine,ivar) = coarseinterp(2,2,ivar)
     &                                + xoff*dudx(i) + yoff*dudy(i)

                 if (ivar.eq.1) then

                    val(ifine,jfine,1)=
     &                val(ifine,jfine,1)-aux(ifine,jfine,1)
		           val(ifine,jfine,1)=dmax1(val(ifine,jfine,1),0.d0)
                 else
                    if (val(ifine,jfine,1).gt.toldry) then
                         velfine(ivar-1)=val(ifine,jfine,ivar)/
     &                                       val(ifine,jfine,1)
                    else
                       redo=.true.
                       velfine(ivar-1)=0.d0
                    endif

                    if ((velfine(ivar-1).ge.coarsevelmax(ivar-1)).or.
     &                (velfine(ivar-1).le.coarsevelmin(ivar-1))) then
                             redo=.true.
                    endif
                 endif
             enddo
             enddo
         enddo

c	   # At this point all of the fine-cells are set.
c	   # If redo=true, then reset the momentum

	   if (redo) then
            do ivar= 2,nvar
                 if (coarseq(2,2,1).gt.toldry) then
                       velcoarse=coarseq(2,2,ivar)/coarseq(2,2,1)
                 else
                       velcoarse=0.d0
                 endif
	           do ico = 1,lratiox
                 do jco = 1,lratioy
                    jfine = (j-2)*lratioy + nghost + jco
                    ifine   = (i-2)*lratiox + nghost + ico
		           val(ifine,jfine,ivar)=velcoarse*val(ifine,jfine,1)
	           enddo
	           enddo
            enddo
	   endif
		


      enddo	!end of coarse loop
      enddo	!end of coarse loop
c
c      if (mcapa .ne. 0) then
c        call fixcapaq(val,aux,mitot,mjtot,valc,auxc,mic,mjc,
c     &                nvar,naux,levc)
c      endif
c
c  overwrite interpolated values with fine grid values, if available.
c
      call intcopy(val,mitot,mjtot,nvar,ilo-nghost,ihi+nghost,
     &             jlo-nghost,jhi+nghost,lev,1,1)

      return
      end
