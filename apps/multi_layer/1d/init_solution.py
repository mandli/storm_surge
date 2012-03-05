#!/usr/bin/env python

import numpy as np

# =============================================================================
#  Bathymetry functions
def set_jump_bathymetry(state,location,depths):
    r"""Set bathymetry such that it has a jump discontinuity
    
    location (float) - X coordinate of jump
    depths (list(float)) - Depth to the left (0) and right (1) of the jump
    """
    x = state.grid.x.centers
    state.aux[0,:] = (x < location) * depths[0] + (x >= location) * depths[1]
    
def set_shelf_bathymetry(state,x0,basin_depth,x1,shelf_depth):
    r"""Set simple shelf like bathymetry
    
    """
    x = state.grid.x.centers
    
    shelf_slope = (basin_depth - shelf_depth) / (x0 - x1)
    
    state.aux[0,:] = (x < x0) * basin_depth
    state.aux[0,:] = (x0 <= x < x1) * (shelf_slope * (x-x0) + basin_depth)
    state.aux[0,:] = (x1 <= x) * shelf_depth
    

def set_h_hat(state,init_location,eta_left,eta_right):
    r"""Set and calculate background states for linearized problem
    
    TODO: Allow this to act over an arbitrary number of layers
    """
     
    x = state.grid.x.centers
    bathy = state.aux[0,:]
    
    left_index = x < init_location
    right_index = x >= init_location
    
    bottom_layer_dry = left_index * (eta_left[1] <= bathy) + right_index * (eta_right[1] <= bathy)
    bottom_layer_wet = left_index * (eta_left[1] > bathy) + right_index * (eta_right[1] > bathy)
    
    state.aux[-2,:] = bottom_layer_wet * (eta_left[0] - eta_left[1])
    state.aux[-1,:] = bottom_layer_wet * (eta_left[1] - bathy)
    state.aux[-2,:] = state.aux[-2,:] + bottom_layer_dry * (eta_left[0] - bathy)
    state.aux[-1,:] = state.aux[-1,:] + bottom_layer_dry * 0.0
    

# =============================================================================
#  Q Init functions
def set_q_quiescent(state):
    r"""Set state variables to a stationary state (quiescent)."""
    
    for layer in xrange(state.problem_data['num_layers']):
        layer_index = 2*layer
        state.q[layer_index,:] = state.aux[layer+3,:] * state.problem_data['rho'][layer]
        state.q[layer_index+1,:] = 0.0
        
    
def set_q_nonzero_velocity(state,location,u_left,u_right):
    r"""Add a jump perturbation to quiescent background state.
    
    init_type == 0
    """
    
    # TODO:  Arbitray number of layer support
    
    # Set quiescent background state
    set_q_quiescent(state)
    
    # Non-zero velocity
    x = state.grid.x.centers
    state.q[1,:] = ((x < init_location) * u_left[0] + (x >= init_location) * u_right[0]) * state.q[0,:]
    state.q[3,:] = ((x < init_location) * u_left[1] + (x >= init_location) * u_right[1]) * state.q[2,:]
    
    
def set_q_simple_wave(state,wave_family,init_location,epsilon):
    r"""Add a perturbation in one wave family to a quiescent background state.
    
    wave_family (int)
    init_location (float)
    epsilon (float)
    
    This used to be init_type == 1
    """
    # TODO:  Arbitray number of layer support
    
    # Set quiescent background state
    set_q_quiescent(state)

    x = state.grid.x.centers
    h_hat = state.aux[-2:,:]
    r = state.problem_data['r']
    g = state.problem_data['g']
    rho = state.problem_data['rho']

    # Riemann problem in one wave family
    gamma = h_hat[1,:] / h_hat[0,:]
    if wave_family == 1:  # Shallow water, left-going
        alpha = 0.5 * (gamma - 1.0 + np.sqrt((gamma - 1.0)**2 + 4.0 * r * gamma))
        eig_value = - np.sqrt(g * h_hat[0,:] * (1.0 + alpha))
    elif wave_family == 2:  # Internal wave, left-going
        alpha = 0.5 * (gamma - 1.0 - np.sqrt((gamma - 1.0)**2 + 4.0 * r * gamma))
        eig_value = - np.sqrt(g * h_hat[0,:] * (1.0 + alpha))
    elif wave_family == 3:  # Internal wave, right-going
        alpha = 0.5 * (gamma - 1.0 - np.sqrt((gamma - 1.0)**2 + 4.0 * r * gamma))
        eig_value = np.sqrt(g * h_hat[0,:] * (1.0 + alpha))
    elif wave_family == 4:  # Shallow water, right-going
        alpha = 0.5 * (gamma - 1.0 + np.sqrt((gamma - 1.0)**2 + 4.0 * r * gamma))
        eig_value = np.sqrt(g * h_hat[0,:] * (1.0 + alpha))

    # Add perturbation
    index = (x < init_location) * (wave_family >= 3)
    state.q[0,index] = state.q[0,index] + rho[0] * epsilon * np.ones((state.q.shape[1]))[index]
    state.q[1,index] = state.q[1,index] + rho[0] * epsilon * eig_value[index]
    state.q[2,index] = state.q[2,index] + rho[1] * epsilon * alpha[index]
    state.q[3,index] = state.q[3,index] + rho[1] * epsilon * alpha[index] * eig_value[index]
    index = (x > init_location) * (wave_family < 3)
    state.q[0,index] = state.q[0,index] + rho[0] * epsilon * np.ones((state.q.shape[1]))[index]
    state.q[1,index] = state.q[1,index] + rho[0] * epsilon * eig_value[index]
    state.q[2,index] = state.q[2,index] + rho[1] * epsilon * alpha[index]
    state.q[3,index] = state.q[3,index] + rho[1] * epsilon * alpha[index] * eig_value[index]
    

def set_q_swe_hump(state,epsilon,location,sigma):
    r"""
    
    init_type == 2
    """
    # TODO:  Arbitray number of layer support
    
    # Set quiescent background state
    set_q_quiescent(state)

    h_hat = state.aux[-2:,:]
    
    gamma = h_hat[1,:] / h_hat[0,:]
    alpha = 0.5 * (gamma - 1.0 + np.sqrt((gamma - 1.0)**2 + 4.0 * r * gamma))
    deta = epsilon * np.exp(-((x - location)/sigma)**2)
    state.q[0,:] = state.q[0,:] + rho[0] * deta
    state.q[2,:] = state.q[2,:] + rho[1] * alpha * deta
            
def set_q_internal_hump(state,epsilon,location,sigma):
    r"""
    
    init_type == 3
    """

    # TODO:  Arbitray number of layer support
    
    # Set quiescent background state
    set_q_quiescent(state)

    h_hat = state.aux[-2:,:]
    
    deta = epsilon * np.exp(-((x - init_location)/sigma)**2)
    state.q[0,:] = state.q[0,:] - rho[0] * deta
    state.q[2,:] = state.q[2,:] + rho[1] * deta
            
def set_q_acta_numerica(state):
    r"""Set q as in the Acta Numerica paper.
    
    init_type == 4
    """

    # TODO:  Arbitray number of layer support
    
    # Set quiescent background state
    set_q_quiescent(state)

    h_hat = state.aux[-2:,:]
    
    gamma = h_hat[1,:] / h_hat[0,:]
    # alpha = 0.5 * (gamma - 1.d0 + np.sqrt((gamma - 1.0)**2+4.0 * r * gamma))
    alpha = np.zeros(())
    xmid = 0.5 * (-180e3 - 80e3)
    deta = (x > -130e3) * (x < -80e3) * epsilon * np.sin((x-xmi)*np.pi / (-80e3 - xmid))
    state.q[2,:] = q[2,:] + rho[1] * alpha * deta
    state.q[0,:] = q[0,:] + rho[0] * deta * (1.0 - alpha)
    

if __name__ == "__main__":
    sys.exit(0)
