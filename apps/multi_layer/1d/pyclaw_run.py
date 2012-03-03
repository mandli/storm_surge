#!/usr/bin/env python
# encoding: utf-8

r""" Run the suite of tests for the 1d two-layer equations"""

import numpy as np

import wind
import init_solution

def before_step(solver,solution,DRY_TOLERANCE=1e-3,RICHARDSON_TOLERANCE=0.95):
    r""""""
    # Extract relevant data
    num_layers = solution.states[0].problem_data['num_layers']
    rho = solution.states[0].problem_data['rho']
    g = solution.states[0].problem_data['g']
    one_minus_r = solution.states[0].problem_data['one_minus_r']
    set_wind = solution.states[0].problem_data['wind_function']
    
    # State arrays
    q = solution.states[0].q
    aux = solution.states[0].aux
    
    # Zero out negative values
    q = q * (q > 0.0)
    
    # Calculate wind
    set_wind(solution)
    
    # Calculate kappa
    h = np.zeros((num_layers,q.shape[1]))
    u = np.zeros(h.shape)
    for layer in xrange(num_layers):
        layer_index = 2*layer
        h[layer,:] = solution.q[layer_index,:] / rho[layer]
        u[layer,:] = (h[layer,:] >= DRY_TOLERANCE) * (q[layer_index+1,:] / q[layer_index,:])
        aux[4,:] = (u[0,:] - u[1,:])**2 / (g * one_minus_r * (h[0,:] + h[1,:]))
    if np.any(aux[4,:] > RICHARDSON_TOLERANCE):
        raise Exception("Richardson tolerance exceeded!")
        
    
def friction_source(solver,state,dt,TOLERANCE=1e-30):
    r""""""
    num_layers = state.problem_data['num_layers']
    manning = state.problem_data['manning']
    g = state.problem_data['g']
    rho = state.problem_data['rho']
    dry_tolerance = state.problem_data['dry_tolerance']
    
    if manning > TOLERANCE:

        for i in xrange(state.q.shape[1]):
            h = state.q[2,i] / rho[1]
            if h < dry_tolerance:
                h = state.q[0,i] / rho[0]
                u = state.q[1,i] / rho[0]
                layer_index = 0
            else:
                u = state.q[2,i] / rho[1]
                layer_index = 1
        
            gamma = u * g * manning**2 / h**(4/3)
            dgamma = 1.0 + dt * gamma
            hu_index = 2 * (layer_index) + 1
            state.q[hu_index,i] = state.q[hu_index,i] / dgamma * rho[layer_index]

def wind_source(solver,state,dt):
    raise NotImplementedError("Wind source terms have not been implemented yet.")


class base_test(object):
    
    def __init__(self):
        
        # Method parameters
        dry_tolerance = 1e-3
        eigen_method = 2
        inundation_method = 2
        entropy_fix = False
        
        # Physical parameters
        manning = 0.025
        num_layers = 2
        
        # Create solver object
        solver = pyclaw.ClawSolver1D()
        solver.bc_lower[0] = 1
        solver.bc_upper[0] = 1
        solver.aux_bc_lower[0] = 1
        solver.aux_bc_upper[0] = 1
        solver.cfl_desired = 0.9
        solver.cfl_max = 1.0
        solver.max_steps = 5000
        
        # Setup Riemann solver
        import riemann
        solver.rp = riemann.rp1_layered_shallow_water
        
        # Set callback functions
        solver.before_step = lambda solver,solution:before_step(solver,solution,
                DRY_TOLERANCE=dry_tolerance,
                RICHARDSON_TOLERANCE=0.95)
        if manning != 0.0:
            solver.step_source = lambda solver,state,dt:friction_source(solver,state,dt,
                TOLERANCE=1e-30)
        
        # General solver parameters (constant for any of the )
        solver.fwave = True
        solver.kernel_language = "Fortran"
        solver.num_waves = num_layers * 2
        solver.limiters = 4
        
        # Create solution
        x = pyclaw.Dimension('x',0.0,1.0,500)
        domain = pyclaw.Domain([x])
        state = pyclaw.State(domain,4,5)

        # Physical parameters
        state.problem_data['num_layers'] = num_layers
        state.problem_data['wind_function'] = wind.set_no_wind
        state.problem_data['g'] = 9.8
        state.problem_data['rho'] = [0.98,1.0]
        state.problem_data['rho_air'] = 1.15e-3
        state.problem_data['manning'] = manning
        state.problem_data['r'] = state.problem_data['rho'][1] / state.problem_data['rho'][0]
        state.problem_data['one_minus_r'] = 1.0 - state.problem_data['r']
        # Solver parameters
        state.problem_data['dry_tolerance'] = dry_tolerance
        state.problem_data['eigen_method'] = eigen_method
        state.problem_data['inundation_method'] = inundation_method
        state.problem_data['entropy_fix'] = entropy_fix
        
        
        # Set simple jump discontinuity in bathy
        init_solution.set_jump_bathymetry(state,0.5,[-1.0,-0.2])
        # Set simple wave initial condition
        init_solution.set_q_simple_wave(state,1,0.45,0.1)
        
        solution = pyclaw.Solution(state,domain)
        solution.t = 0.0
        
        self.claw = pyclaw.Controller()
        self.claw.tfinal = 1.0
        self.claw.solution = solution
        self.claw.solver = solver
        self.claw.num_output_times = 10
        self.claw.outdir = "_output"
            
    def run(self):
        return self.claw.run()
        
        
if __name__ == "__main__":
    use_petsc = False
    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw
        
    test = base_test()
    test.run()