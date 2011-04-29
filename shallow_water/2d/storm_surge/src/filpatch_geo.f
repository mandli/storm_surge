c
c ---------------------------------------------------------------
c
        recursive subroutine filrecur(level,nvar,valbig,aux,naux,
     1                      time,mitot,mjtot,
     2                      nrowst,ncolst,ilo,ihi,jlo,jhi)

c :::::::::::::::::::::::::::: FILPATCH :::::::::::::::::::::::::;
c
c  fill the portion of valbig from rows  nrowst
c                             and  cols  ncolst
c  the patch can also be described by the corners (xlp,ybp) by (xrp,ytp).
c  vals are needed at time time, and level level,
c
c  first fill with  values obtainable from the level level
c  grids. if any left unfilled, then enlarge remaining rectangle of
c  unfilled values by 1 (for later linear interp), and recusively
c  obtain the remaining values from  coarser levels.
c
c :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;
    
      use geoclaw_module
    
      implicit double precision (a-h,o-z)

      include  "call.i"

      logical   set, sticksout
      dimension valbig(mitot,mjtot,nvar), aux(mitot,mjtot,naux)

      dimension coarseq(3,3,3)
      dimension coarseeta(3,3)
      dimension coarseb(3,3)
      dimension coarseinterp(3,3,3)
      dimension fineval(3)
      dimension velfine(2)
      dimension coarsevelmax(2)
      dimension coarsevelmin(2)

      logical redo

c  use stack-based scratch arrays instead of alloc, since dont really
c  need to save beyond these routines, and to allow dynamic memory resizing
c
c     use 1d scratch arrays that are potentially the same size as 
c     current grid, since may not coarsen.
c     need to make it 1d instead of 2 and do own indexing, since
c     when pass it in to subroutines they treat it as having different
c     dimensions than the max size need to allocate here
c
      dimension valcrse((ihi-ilo+2)*(jhi-jlo+2)*nvar)  ! NB this is a 1D array 
      dimension auxcrse((ihi-ilo+2)*(jhi-jlo+2)*naux)  ! the +2 is to expand on coarse grid to enclose fine
c
      dimension flaguse(ihi-ilo+1,jhi-jlo+1)


c      iadflag(i,j)    =  locuse + i-1+(j-1)*nrowp
c      ivalc(i,j,ivar) =  loccrse + (i - 1) + nrowc*(j - 1)
c     &                     + nrowc*ncolc*(ivar-1)
      ivalc(i,j,ivar) = i + nrowc*(j - 1)
     &                    + nrowc*ncolc*(ivar-1)
c
c     # index into first component of aux = topo:
c      iauxc(i,j) =  locauxc + (i - 1) + nrowc*(j - 1)
      iauxc(i,j) =  i + nrowc*(j - 1)

      sticksout(iplo,iphi,jplo,jphi)  =
     &            (iplo .lt. 0 .or. jplo .lt. 0 .or.
     &             iphi .ge. iregsz(levc) .or. jphi .ge. jregsz(levc))


!--      write(*,*)" entering filrecur with level ",level
!--      write(*,*)"     and patch indices ilo,ihi,jlo,jhi ",
!--     &             ilo,ihi,jlo,jhi

c         write(*,*)" in filrecur for level ",level,mitot,mjtot
c
c We begin by filling values for grids at level level. If all values can be
c filled in this way, we return;

        nrowp   = ihi - ilo + 1
        ncolp   = jhi - jlo + 1
c        locuse  = igetsp(nrowp*ncolp)
        hxf     = hxposs(level)
        hyf     = hyposs(level)
        xlp     = xlower + ilo*hxf
        xrp     = xlower + (ihi+1)*hxf
        ybp     = ylower + jlo*hyf
        ytp     = ylower + (jhi+1)*hyf

        call intfil
     &  (valbig,mitot,mjtot,time,flaguse,nrowst,ncolst,
     &   ilo,ihi,jlo,jhi,level,nvar,naux)
c     &  (valbig,mitot,mjtot,time,locuse,nrowst,ncolst,

c
c Trimbd returns set = true if all of the entries are filled (=1.).
c set = false, otherwise. If set = true, then no other levels are
c are required to interpolate, and we return.
c
c Note that the used array is filled entirely in intfil, i.e. the
c marking done there also takes into account the points filled by
c the boundary conditions. bc2amr will be called later, after all 4
c boundary pieces filled.

c        call trimbd(alloc(locuse),nrowp,ncolp,set,il,ir,jb,jt)
        call trimbd(flaguse,nrowp,ncolp,set,il,ir,jb,jt)

        if (set) go to 90 ! all done except for bcs
c
c otherwise make recursive calls to coarser levels to fill remaining unset points
c
        if (level .eq. 1) then
           write(outunit,*)" error in filrecur - level 1 not set"
           write(outunit,900) nrowst,ncolst
           write(*,*)" error in filrecur - level 1 not set"
           write(*,*)" should not need more recursion "
           write(*,*)" to set patch boundaries"
           write(*,900) nrowst,ncolst
900        format("start at row: ",i4," col ",i4)
           stop
        endif

c set = false. we will have to interpolate some values from coarser
c levels. We begin by initializing the level level arrays, so that we can use
c purely recursive formulation for interpolating.


        levc = level - 1
        hxc  = hxposs(levc)
        hyc  = hyposs(levc)

        isl  = il + ilo - 1
        isr  = ir + ilo - 1
        jsb  = jb + jlo - 1
        jst  = jt + jlo - 1
c
c       coarsen
        lratiox = intratx(levc)
        lratioy = intraty(levc)
        iplo   = (isl-lratiox  +nghost*lratiox)/lratiox - nghost
        jplo   = (jsb-lratioy  +nghost*lratioy)/lratioy - nghost
        iphi   = (isr+lratiox  )/lratiox
        jphi   = (jst+lratioy  )/lratioy

        xlc  =  xlower + iplo*hxc
        ybc  =  ylower + jplo*hyc
        xrc  =  xlower + (iphi+1)*hxc
        ytc  =  ylower + (jphi+1)*hyc

        nrowc   =  iphi - iplo + 1
        ncolc   =  jphi - jplo + 1
        ntot    = nrowc*ncolc*(nvar+naux)
c        write(*,876) nrowc,ncolc, ihi-ilo+2,jhi-jlo+2
 876    format(" needed coarse grid size ",2i5," allocated ",2i5)
        if (nrowc .gt. ihi-ilo+2 .or. ncolc .gt. jhi-jlo+2) then
            write(*,*)" did not make big enough work space in filrecur "
            write(*,*)" need coarse space with nrowc,ncolc ",nrowc,ncolc
            write(6,*)" made space for ilo,ihi,jlo,jhi ",ilo,ihi,jlo,jhi
            stop
        endif
c        loccrse = igetsp(ntot)
c        locauxc = loccrse + nrowc*ncolc*nvar
        if (naux.gt.0) then
              maxmx = nrowc - 2*nghost
              mx = maxmx
              maxmy = ncolc - 2*nghost
              my = maxmy
              xl = xlc + nghost*hxc
              yb = ybc + nghost*hyc
              call setaux(maxmx,maxmy,nghost,mx,my,xl,yb,hxc,hyc,
     &                    naux,auxcrse)
c     &                    naux,alloc(locauxc))
        endif

        if ((xperdom .or. (yperdom .or. spheredom)) .and.
     &       sticksout(iplo,iphi,jplo,jphi)) then
            call prefilrecur(levc,nvar,valcrse,auxcrse,
     1                    naux,time,nrowc,ncolc,1,1,
     2                    iplo,iphi,jplo,jphi)
        else
c          call filpatch2(levc,nvar,alloc(loccrse),alloc(locauxc),naux,
          call filrecur(levc,nvar,valcrse,auxcrse,naux,
     1                   time,nrowc,ncolc,1,1,
     2                   iplo,iphi,jplo,jphi)
        endif

c       interpolate back up

20      continue


c       #loop through all fine cells interpolating eta,hu,hv

        toldry= drytolerance

        do 100 iff = 1,nrowp
          ic = 2 + (iff - (isl - ilo) - 1)/lratiox
          eta1 = (-0.5d0+dble(mod(iff-1,lratiox)))/dble(lratiox)

        do 100 jf  = 1,ncolp
          jc = 2 + (jf -(jsb-jlo)-1)/lratioy
          eta2 = (-0.5d0+dble(mod(jf -1,lratioy)))/dble(lratioy)

c           flag = alloc(iadflag(iff,jf))
           flag = flaguse(iff,jf)
           if (flag .eq. 0.0) then


c          # 5 coarse values surrounding and including the current coarse value are set
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
                  coarseb(2+ii,2+jj)= auxcrse(iauxc(ic+ii,jc+jj))
                  do ivar=1,nvar
                   coarseq(2+ii,2+jj,ivar)=
     &                            valcrse(ivalc(ic+ii,jc+jj,ivar))
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
     &                 coarseinterp(3,2,ivar)-coarseinterp(1,2,ivar))
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

                   slopex=slopex*sxlimit
                   slopey=slopey*sylimit

                   fineval(ivar)=coarseinterp(2,2,ivar) + eta1*slopex
     &                  + eta2*slopey
                   if (ivar.eq.1) then

c               ### if coarse cell is dry then force fine grid to use
c               ### coarse bathymetry. avoids inconsistency
                      if (coarseq(2,2,1) .lt. toldry) then
c                        ### debug output:
c                        write(6,*)"setting to coarse bathy for ",
c    .                        iff+nrowst-1,jf+ncolst-1
c                        write(6,*)"    old bathy ",
c    .                           aux(iff+nrowst-1,jf+ncolst-1,1)
c                        write(6,*)"    new bathy ",coarseb(2,2)

                         aux(iff+nrowst-1,jf+ncolst-1,1)=coarseb(2,2)
                      endif

                      fineval(1)=fineval(1)
     &                     -aux(iff+nrowst-1,jf+ncolst-1,1)
                      fineval(1)=max(fineval(1),0.d0)
                   else
                      redo=.false.
                      if (fineval(1).gt.toldry) then
                         velfine(ivar-1)=fineval(ivar)/fineval(1)
                      else
                         redo=.true.
                         velfine(ivar-1)=0.d0
                      endif

                      if ((velfine(ivar-1).gt.coarsevelmax(ivar-1).or.
     &                  velfine(ivar-1).lt.coarsevelmin(ivar-1))
     &                            .or.redo) then
                          if (coarseq(2,2,1).gt.toldry) then
                           fineval(ivar)=fineval(1)*coarseq(2,2,ivar)/
     &                         coarseq(2,2,1)
                          else
                           fineval(ivar)=0.d0

                          endif
                      endif
                   endif
                  valbig(iff+nrowst-1,jf+ncolst-1,ivar)=fineval(ivar)
                enddo

         endif  ! end if flag
 100     continue

c        call reclam(loccrse,ntot)

 90       continue
c
c  set bcs, whether or not recursive calls needed. set any part of patch that stuck out
c
        call bc2amr(valbig,aux,mitot,mjtot,nvar,naux,
     1              hxf,hyf,level,time,
     2              xlp,xrp,ybp,ytp,
     3              xlower,ylower,xupper,yupper,
     4              xperdom,yperdom,spheredom)


c        call reclam(locuse,nrowp*ncolp)

        return
        end
