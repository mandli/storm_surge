! ============================================================================
!  Program:     /Users/mandli/src/local/shallow_water/2d/src
!  File:        riemann_fwave_stormsurge
!  Created:     2010-06-18
!  Author:      Kyle Mandli
! ============================================================================
!      Copyright (C) 2010-06-18 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================
!  This is a modification of the riemann solver riemann_fwave from the 
!  geoclaw lib (ver 4.5.0).  It includes the wind source term for storm 
!  surges.  (KTM 2010-06-18)
! ============================================================================

subroutine riemann_fwave_stormsurge(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,s1,s2,drytol,tau,wind,g,sw,fw)

      ! solve shallow water equations given single left and right states
      ! solution has two waves.
      ! flux - source is decomposed.

      use hurricane_module

      implicit none

      !input
      integer, intent(in) :: meqn,mwaves

      double precision, intent(in) :: hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,s1,s2
      double precision, intent(in) :: hvL,hvR,vL,vR
      double precision, intent(in) :: drytol,tau,wind,g

      double precision, intent(out) :: sw(mwaves)
      double precision, intent(out) :: fw(meqn,mwaves)

      !local
      double precision :: delh,delhu,delphi,delb,delhdecomp,delphidecomp
      double precision :: deldelh,deldelphi
      double precision :: beta1,beta2
      
      integer icom,jcom
      double precision dtcom,dxcom,dycom,tcom
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      sw = 0.d0
      fw = 0.d0
      
      return

      !determine del vectors
      delh = hR-hL
      delhu = huR-huL
      delphi = phiR-phiL-dxcom*tau*wind
      delb = bR-bL

      deldelphi = -g*0.5d0*(hR+hL)*delb
      delphidecomp = delphi - deldelphi

      !flux decomposition
      beta1 = (s2*delhu - delphidecomp)/(s2-s1)
      beta2 = (delphidecomp - s1*delhu)/(s2-s1)

      sw(1)=s1
      sw(2)=0.5d0*(s1+s2)
      sw(3)=s2
      ! 1st nonlinear wave
      fw(1,1) = beta1
      fw(2,1) = beta1*s1
      fw(3,1) = beta1*vL
      ! 2nd nonlinear wave
      fw(1,3) = beta2
      fw(2,3) = beta2*s2
      fw(3,3) = beta2*vR
      ! advection of transverse wave
      fw(1,2) = 0.d0
      fw(2,2) = 0.d0
      fw(3,2) = hR*uR*vR - hL*uL*vL -fw(3,1)-fw(3,3)
      return

      end