#!/usr/bin/env python
# encoding: utf-8

r""" Run the suite of tests for the 1d two-layer equations"""

import os
import sys
import types

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
    # h = np.zeros((num_layers,q.shape[1]))
    # u = np.zeros(h.shape)
    # for layer in xrange(num_layers):
    #     layer_index = 2*layer
    #     h[layer,:] = solution.q[layer_index,:] / rho[layer]
    #     u[layer,:] = (h[layer,:] >= DRY_TOLERANCE) * (q[layer_index+1,:] / q[layer_index,:])
    #     aux[4,:] = (u[0,:] - u[1,:])**2 / (g * one_minus_r * (h[0,:] + h[1,:]))
    # if np.any(aux[4,:] > RICHARDSON_TOLERANCE):
    #     raise Exception("Richardson tolerance exceeded!")
        
    
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


def ml_swe(use_petsc=False,iplot=False,htmlplot=False,outdir='./_output',solver_type='classic'):
    r""""""
    
    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw
        
    # Method parameters
    dry_tolerance = 1e-3
    eigen_method = 2
    inundation_method = 2
    entropy_fix = False
        
    # Physical parameters
    manning = 0.025
    num_layers = 2
    rho = [0.95,1.0]
    
    # Init condition
    wave_family = 4
        
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
    state.problem_data['rho'] = rho
    state.problem_data['rho_air'] = 1.15e-3
    state.problem_data['manning'] = manning
    state.problem_data['r'] = state.problem_data['rho'][0] / state.problem_data['rho'][1]
    state.problem_data['one_minus_r'] = 1.0 - state.problem_data['r']

    # Solver parameters (needed for the Riemann solver)
    state.problem_data['dry_tolerance'] = dry_tolerance
    state.problem_data['eigen_method'] = eigen_method
    state.problem_data['inundation_method'] = inundation_method
    state.problem_data['entropy_fix'] = entropy_fix
        
    # Set simple jump discontinuity in bathy
    bathy_type = 1
    bathy_location = 0.5
    init_solution.set_jump_bathymetry(state,bathy_location,[-1.0,-0.2])
        
    # Set background states
    init_location = 0.45
    eta_left = [0.0,-0.6]
    eta_right = [0.0,-0.6]
    init_solution.set_h_hat(state,init_location,eta_left,eta_right)
    
    # Set perturbation in particuar wave family
    # init_solution.set_q_quiescent(state)
    init_solution.set_q_simple_wave(state,wave_family,init_location,0.1)
        
    solution = pyclaw.Solution(state,domain)
    solution.t = 0.0
        
    # Create controller
    claw = pyclaw.Controller()
    claw.solution = solution
    claw.solver = solver
    # claw.output_style = 1
    # claw.tfinal = 0.5
    claw.output_style = 3
    claw.nstepout = 1
    claw.num_output_times = 50
    claw.outdir = outdir
    claw.write_aux_init = True

    state = claw.run()
    
    if htmlplot or iplot:
        import pyclaw.plot
    
        setplot_path = './setplot.py'
        plot_kargs = {'bathy_type':bathy_type,'bathy_location':bathy_location,
                      'init_type':1,'wave_family':wave_family,
                      'dry_tolerance':dry_tolerance,
                      'rho':rho}
                      #'x0','x1',
    
        # Grab and import the setplot function
        path = os.path.abspath(os.path.expandvars(os.path.expanduser(setplot_path)))
        setplot_module_dir = os.path.dirname(path)
        setplot_module_name = os.path.splitext(os.path.basename(setplot_path))[0]
        sys.path.insert(0,setplot_module_dir)
        setplot_module = __import__(setplot_module_name)
        reload(setplot_module)
        setplot = lambda plotdata:setplot_module.setplot(plotdata,**plot_kargs)
        
        if not isinstance(setplot,types.FunctionType):
            raise ImportError("Failed importing %s.setplot" % setplot_module_name)
        

        if iplot:     
            from visclaw import Iplotclaw
        
            ip=Iplotclaw.Iplotclaw(setplot=setplot)
            ip.plotdata.outdir = outdir
            ip.plotdata.format = 'ascii'
        
            ip.plotloop()
            
        if htmlplot:  
            from visclaw import plotclaw
            
            plotclaw.plotclaw(outdir,format='ascii',setplot=setplot)
        
    return claw.solution.q
        
        
if __name__ == "__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(ml_swe)