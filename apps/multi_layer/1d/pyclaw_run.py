#!/usr/bin/env python
# encoding: utf-8

r""" Run the suite of tests for the 1d two-layer equations"""

def set_wind(solution):
    
    if wind_type == 0:
        solution.aux[1,:] = 0.0
    if wind_type == 1:

def before_step(solver,solution):
    
    # Extract relevant data
    num_layers = solution.states[0].prob_data['num_layers']
    rho = solution.states[0].prob_data['rho']
    g = solution.states[0].prob_data['g']
    
    # Zero out negative values
    solution.q = solution.q * (solution.q > 0.0)
    
    # Calculate wind
    set_wind(solution)
    
    # Calculate kappa
    for layer in xrange(num_layers):
        layer_index = 2*layer
        h[layer,:] = solution.q[layer_index,:] / rho[layer]
        u = (h >= DRY_TOLERANCE) * (solution.q[layer_index+1,:] / solution.q[layer_index,:])
        aux[4,:] = (u[0] - u[1])**2 / (g * one_minus_r * (h[0,:] + h[1,:]))
    if np.any(aux[4,:] > RICHARDSON_TOLERANCE):
        raise Exception("Richardson tolerance exceeded!")

class base_test(object):
    
    def __init__(self,mx=500,eigen_method=4,inundation_method=2,layers=2):
        
        # Create solver object
        solver = pyclaw.Claw()
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
        rp_data = solver.rp.rp1_layer_shallow_water_module
        rp_data.dry_tolerance = 1e-3
        rp_data.eigen_method = 4
        rp_data.inundation_method = 1
        rp_data.entropy_fix = False
        
        # General solver parameters (constant for any of the )
        solver.solver.fwave = True
        solver.kernel_language = "Fortran"
        solver.num_waves = self.layers * 2
        solver.limiters = 4
        
        # Algorithm parameters
        eigen_method = 4
        inundation_method = 1
        entropy_fix = False
        richardson_tolerance = 0.95
        wave_tolerance = [1e-1,2e-1] # AMR parameter
        dry_tolerance = 1e-3
        dry_limit = False
        
        # Create solution
        x = pyclaw.Dimension('x',0.0,1.0,mx)
        domain = pyclaw.Domain([x])
        state = pyclaw.State(domain,4,5)

        # Physical parameters
        state.problem_data['g'] = 9.8
        state.problem_data['rho'] = [1025.,1028.]
        state.problem_data['rho_air'] = 1.15e-3
        state.problem_data['manning'] = 0.025
        state.problem_data['r'] = state.problem_data['rho'][1] / state.problem_data['rho'][0]
        state.problem_data['one_minus_r'] = 1.0 - state.problem_data['r']
        
        solution = pyclaw.Solution(domain,state)
        solution.t = 0.0
        init_solution(solution)
        
        claw = pyclaw.Controller()
        claw.tfinal = 1.0
        claw.solution = solution
        claw.solver = solver
        claw.num_output_times = 10
        claw.outdir = ""
            
    def run(self):
        return claw.run()
        
def run_multilayer_swe(mx=500,r=0.98,use_petsc=False):
    r""""""
    import numpy as np
    
    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw
        
if __name__ == "__main__":
    pass