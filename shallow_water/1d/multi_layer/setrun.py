""" 
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

import os
import numpy as np

import pyclaw.data as data
from multilayer_data import MultilayerData

wave_family = 3

#------------------------------
def setrun(claw_pkg='Classic'):
#------------------------------
    
    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "Classic4" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """ 
    
    assert claw_pkg.lower() == 'classic',  "Expected claw_pkg = 'classic'"

    rundata = data.ClawRunData(pkg=claw_pkg, ndim=1)
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.ndim = 1
    
    # Lower and upper edge of computational domain:
    clawdata.xlower = 0.0
    clawdata.xupper = 1.0
        
    # Number of grid cells:
    clawdata.mx = 500

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 4

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 5
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.mcapa = 0
    
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.outstyle = 1
    
    # For outstyle 1 specified in hours
    num_hours = 40
    step = 0.25

    if clawdata.outstyle==1:
        # Output nout frames at equally spaced times up to tfinal:
        # clawdata.nout = int(num_hours / step) + int(np.ceil(RAMP_UP_TIME / (step*60**2)))
        # clawdata.tfinal = num_hours * 60.0**2
        # clawdata.nout = 40
        # clawdata.tfinal = 40
        clawdata.nout = 40
        if wave_family < 4 and wave_family > 1:
            clawdata.tfinal = 1.0
        else:
            clawdata.tfinal = 0.1
        

    elif clawdata.outstyle == 2:
        # Specify a list of output times.  
        clawdata.tout =  [0.5, 1.0]   # used if outstyle == 2
        clawdata.nout = len(clawdata.tout)

    elif clawdata.outstyle == 3:
        # Output every iout timesteps with a total of ntot time steps:
        iout = 1
        ntot = 200
        clawdata.iout = [iout, ntot]
    


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = 1
    
    # Initial time step for variable dt.  
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.575e-03
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used, and max to allow without 
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.0
    
    # Maximum number of time steps to allow between output times:
    clawdata.max_steps = 5000

    
    

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Transverse order for 2d or 3d (not used in 1d):
    clawdata.order_trans = 0
    
    # Number of waves in the Riemann solution:
    clawdata.mwaves = 4
    
    # List of limiters to use for each wave family:  
    # Required:  len(mthlim) == mwaves
    clawdata.mthlim = [3,3,3,3]
    
    # Source terms splitting:
    #   src_split == 0  => no source term (src routine never called)
    #   src_split == 1  => Godunov (1st order) splitting used, 
    #   src_split == 2  => Strang (2nd order) splitting used,  not recommended.
    clawdata.src_split = 0
    
    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.mbc = 2
    
    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity
    
    clawdata.mthbc_xlower = 1
    clawdata.mthbc_xupper = 1
    
    return rundata
    # end of function setrun
    # ----------------------

def set_multilayer_data():
    
    prob_data = MultilayerData()
    
    # Physics
    prob_data.rho_air = 1.15e-3
    prob_data.rho_1 = 0.98
    prob_data.rho_2 = 1.0
    
    # Algorithm
    prob_data.dry_tolerance = 1e-3
    prob_data.eigen_method = 4
    
    # Initial condition
    prob_data.init_type = 1
    prob_data.init_location = 0.45
    prob_data.wave_family = wave_family
    prob_data.eta_1 = 0.0
    prob_data.eta_2 = -0.6
    prob_data.epsilon = 0.1
    prob_data.sigma = 0.02
    
    # Bathymetry
    prob_data.bathy_location = 0.5
    prob_data.bathy_left = -1.0
    prob_data.bathy_right = -0.2
    
    # Wind
    prob_data.wind_type = 1
    prob_data.A = 5.0
    prob_data.omega = 2.0
    prob_data.N = 2.0
    prob_data.t_length = 10.0
    prob_data.B = 1.5
    prob_data.Pn = 1005.0
    prob_data.Pc = 950.0
    prob_data.hurricane_velocity = 5.0
    prob_data.R_eye_init = 0.0
    prob_data.ramp_up_time = 0.0
    
    return prob_data

if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
    	rundata = setrun()
    prob_data = set_multilayer()

    rundata.write()
    prob_data.write()
