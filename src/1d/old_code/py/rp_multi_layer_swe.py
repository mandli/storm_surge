#!/usr/bin/env python
# encoding: utf-8
r"""Reimann solvers to solve the multilayer shallow water system

..math:
    
    (h_1)_t + (h_1 u_1)_x = 0
    (h_1 u_1)_t + (h_1 u_1^2 + 1/2 g h_1^2)_x = -gh_1(h_2)_x - gh_1b_x
    (h_2)_t + (h_2 u_2)_x = 0
    (h_2 u_2)_t + (h_2 u_2^2 + 1/2 g h_2^2)_x = -gh_2(h_1)_x - gh_2b_x
    
All solvers expect *aux_global* to contain
 - *g* - (float) Gravitational acceleration constant
 - *r* - (float) Ratio of the densities of the layers

:Authors:
    Kyle T. Mandli (2009-9-2):  Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import numpy as np
import scipy.linalg as linalg

def rp_multi_layer_swe_1d_linearized(q_l,q_r,aux_l,aux_r,aux_global):
    r"""Multilayer swe solver in 1d using a linearized f-wave approach
    
    Solve the shallow water equations using the non-conservative system with
    a flux matrix incorporating both the bathymetry term and the 
    non-conservative products.  See the notes found in this directory for more
    explanation.
    
    :Version: 0.1 (2010-9-21)
    """
    
    # Solver constants
    meqn = 4
    mwaves = meqn
    DRY_TOLERANCE = 1e-3
    
    # Arrays we are returning
    nrp = np.size(q_l,0)
    fwave = np.empty( (nrp, meqn, mwaves) )
    s = np.empty( (nrp, mwaves) )
    amdq = np.empty( (nrp, meqn) )
    apdq = np.empty( (nrp, meqn) )
    
    # Work arrays
    delta = np.empty(4)
    R = np.empty((4,4))
    beta = np.empty(4)
    
    for i in xrange(nrp):
        # Perform some pre-processing on the state variables
        
        h1_l = q_l[i,0]
        h1_r = q_r[i,0]
        if h1_l > DRY_TOLERANCE:
            u1_l = q_l[i,1] / q_l[i,0]
        else:
            h1_l = 0.0
            u1_l = 0.0
        if h1_r > DRY_TOLERANCE:
            u1_r = q_r[i,1] / q_l[i,0]
        else:
            h1_r = 0.0
            u1_r = 0.0
        
        h2_l = q_l[i,3]
        h2_r = q_r[i,3]
        if h2_l > DRY_TOLERANCE:
            u2_l = q_l[i,4] / q_l[i,3]
        else:
            h2_l = 0.0
            u2_l = 0.0
        if h2_r > DRY_TOLERANCE:
            u2_r = q_r[i,4] / q_l[i,3]
        else:
            h2_r = 0.0
            u2_r = 0.0
                    
        # Calculate the speeds, out to order 1-r
        # External speed, left
        s[i,0] = -np.sqrt(g*(q_l[i,0]+q_l[i,3])) + q_l[i,0] * q_l[i,3] / (2.0 * (q_l[i,0]+q_l[i,3])**(3/2)) * (1.0-r)
        # Internal speed, left
        s[i,1] = -np.sqrt(g*(q_l[i,0] * q_l[i,3]) / (q_l[i,0] + q_l[i,3]) * (1.0-r))
        # Internal speed, right
        s[i,2] = np.sqrt(g*(q_r[i,0] * q_r[i,3]) / (q_r[i,0] + q_r[i,3]) * (1.0-r))
        # External speed, right
        s[i,3] = np.sqrt(g*(q_r[i,0]+q_r[i,3])) - q_r[i,0] * q_r[i,3] / (2.0 * (q_r[i,0]+q_r[i,3])**(3/2)) * (1.0-r)
                    
        # Calculate the delta vector
        delta = 0.0
        delta[0] = q_r[i,1] - q_l[i,1]
        if q_r[i,0] > DRY_TOLERANCE:
            delta[1] = delta[1] + q_r[i,1]**2 / q_r[i,0]
        delta[1] = delta[1] + 0.5 * g * q_l[i,0]**2
        if q_l[i,0] > DRY_TOLERANCE:
            delta[1] = delta[1] - q_l[i,1]**2 / q_l[i,0] 
        delta[1] = delta[1] - 0.5 * g * q_l[i,0]**2
        delta[2] = q_r[i,4] - q_l[i,4]
        if q_r[i,3] > DRY_TOLERANCE:
            delta[4] = delta[4] + q_r[i,4]**2 / q_r[i,3]
        delta[4] = delta[4] + 0.5 * g * q_l[i,3]**2
        if q_l[i,3] > DRY_TOLERANCE:
            delta[4] = delta[4] - q_l[i,4]**2 / q_l[i,3]
        delta[4] = delta[4] - q_r[i,4]**2 / q_r[i,3]
        
        # Calculate the Eigen-vector matrix
        R = 0.0
        R[0,:] = 1.0
        R[1,:] = s[i,:]
        R[2,:] = [(s[i,0]**2 - g*q_l[i,0]) / (g*r*q_l[i,0]),
                  (s[i,1]**2 - g*q_l[i,0]) / (g*r*q_l[i,0]),
                  (s[i,2]**2 - g*q_r[i,0]) / (g*r*q_r[i,0]),
                  (s[i,3]**2 - g*q_r[i,0]) / (g*r*q_r[i,0])]
        R[3,:] = R[2,:] * s[i,:]
    
        # Solve the system to find the wave strengths
        try:
            beta = linalg.solve(R,delta)
        except LinAlgError:
            print "ERROR:  Linear solver found singular system."
            print "R = %s" % R
            print "delta = %s" % delta
            raise
        
        # Compute each wave
        fwave[i,:,0] = beta[0] * R[:,0]
        fwave[i,:,1] = beta[1] * R[:,1]
        fwave[i,:,2] = beta[2] * R[:,2]
        fwave[i,:,3] = beta[3] * R[:,3]
        
        # Compute fluxes
        for m in xrange(meqn):
            for mw in xrange(mwaves):
                amdq[i,m] += np.min(s[i,mw],0.0) * fwave[i,m,mw]
                apdq[i,m] += np.max(s[i,mw],0.0) * fwave[i,m,mw]
    
    return fwave,s,amdq,apdq

def rp_multi_layer_swe_1d_relax(q_l,q_r,aux_l,aux_r,aux_global):
    r"""Multilayer shallow water solver in 1d using a relaxation system
    
    Riemann solver for the relaxation system described in

    Abgrall and Karni. Two-Layer Shallow Water Systems: A Relaxation Approach.
    SIAM Journal of Scientific Computing (2009) vol. 31 (3) pp. 1603-1627

    which solves

    .. math:

        (\rho_1 h_1)_t + (\rho_1 h_1 u_1)_x = 0
        (\rho_1 h_1 u_1)_t + (\rho_1 h_1 u_1^2 + 1/2 g \rho_1 h_1^2 + \rho_2 h_1 h_2')_x = \rho_2 g h_2 (h_1')_x - g \rho_1 h_1 b_x
        (\rho_2 h_2)_t + (\rho_2 h_2 u_2)_x = 0
        (\rho_2 h_2 u_2)_t + (\rho_2 h_2 u_2^2 + 1/2 g \rho_2 h_2^2)_x = -\rho_2 g h_2 (h_1')_x - g \rho_2 h_2 b_x
        (h_1')_t + U_1^* (h_1')_x = \frac{h_1' - h_1}{\epsilon}
        (h_2')_t + U_2^* (h_2')_x = \frac{h_2' - h_2}{\epsilon}
        
    In this solver, we ignore the relaxation terms on the right hand side of 
    equation 5 and 6.
     
    :Version: 0.1 (2009-9-2)
    """
    
    # Solver constants
    meqn = 6
    mwaves = meqn
    
    # Arrays we are returning
    nrp = np.size(q_l,0)
    wave = np.empty( (nrp, meqn, mwaves) )
    s = np.empty( (nrp, mwaves) )
    amdq = np.empty( (nrp, meqn) )
    apdq = np.empty( (nrp, meqn) )
    
    c1 = g*h1+r*g*h2_prime
    c2 = g*h2
    U1_star = u1
    U2_star = u2
    
    s[:,0] = u1-c1
    s[:,1] = u1+c1
    s[:,2] = u2-c2
    s[:,3] = u2+c2
    s[:,4] = U1_star
    s[:,5] = U2_star
    
    return wave, s, amdq, apdq