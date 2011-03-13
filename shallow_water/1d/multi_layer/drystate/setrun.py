""" 
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

import os
import numpy as np

import pyclaw.data as data

# Simple hurricane data format
class MultilayerData(data.Data):
    def __init__(self,rundata):
        super(MultilayerData,self).__init__()

        # # Physics
        # self.add_attribute('r',0.98)
        # 
        # # Algorithm
        # self.add_attribute('dry_tolerance',1e-3)
        # 
        # # Initial condition
        # self.add_attribute('wave_family',4)
        # self.add_attribute('eta_1',0.0)
        # self.add_attribute('eta_2',-200.0)
        # self.add_attribute('epsilon',1.0)
        # self.add_attribute('init_location',500.0)
        # 
        # # Bathymetry
        # self.add_attribute('bathy_location',1000.0)
        # self.add_attribute('bathy_left',-2000.0)
        # self.add_attribute('bathy_right',-100.0)

        # Physics
        self.add_attribute('rho_1',0.95)
        self.add_attribute('rho_2',1.0)
        
        # Algorithm
        self.add_attribute('dry_tolerance',1e-3)
        
        # Initial condition
        self.add_attribute('wave_family',4)
        self.add_attribute('eta_1',0.0)
        # self.add_attribute('eta_2',-1.0)
        # self.add_attribute('eta_2',-0.8)
        self.add_attribute('eta_2',-0.6)
        self.add_attribute('epsilon',0.05)
        self.add_attribute('init_location',0.48)
        
        # Bathymetry
        self.add_attribute('bathy_location',0.5)
        self.add_attribute('bathy_left',-1.0)
        self.add_attribute('bathy_right',-1.0)
        # self.add_attribute('bathy_right',-0.2)
        
        
    def write(self,file='./problem.data',datasource='setrun.py'):
        """Write out the data file to the path given"""

        print "Creating data file %s" % file
        out_file = data.open_datafile(file)
        
        data.data_write(out_file,self,'rho_1','(Density of top layer)')
        data.data_write(out_file,self,'rho_2','(Density of bottom layer)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'dry_tolerance','(Dry state tolerance)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'init_location','(Location for perturbation)')
        data.data_write(out_file,self,'wave_family','(Wave family of the perturbation)')
        data.data_write(out_file,self,'eta_1','(Steady state top surface)')
        data.data_write(out_file,self,'eta_2','(Steady state internal surface)')
        data.data_write(out_file,self,'epsilon','(Perturbation strength)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'bathy_location','(Bathymetry jump location)')
        data.data_write(out_file,self,'bathy_left','(Depth to left of bathy_location)')
        data.data_write(out_file,self,'bathy_right','(Depth to right of bathy_location)')
        
        out_file.close()



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
    # clawdata.xlower = 0.0
    # clawdata.xupper = 2.0e3
        
    # Number of grid cells:
    clawdata.mx = 500

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 4

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 2
    
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

    clawdata.outstyle = 3
    
    # For outstyle 1 specified in hours
    num_hours = 40
    step = 0.25

    if clawdata.outstyle==1:
        # Output nout frames at equally spaced times up to tfinal:
        # clawdata.nout = int(num_hours / step) + int(np.ceil(RAMP_UP_TIME / (step*60**2)))
        # clawdata.tfinal = num_hours * 60.0**2
        # clawdata.nout = 40
        # clawdata.tfinal = 40
        clawdata.nout = 20
        clawdata.tfinal = 0.5

    elif clawdata.outstyle == 2:
        # Specify a list of output times.  
        clawdata.tout =  [0.5, 1.0]   # used if outstyle == 2
        clawdata.nout = len(clawdata.tout)

    elif clawdata.outstyle == 3:
        # Output every iout timesteps with a total of ntot time steps:
        iout = 1
        ntot = 40
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
    clawdata.max_steps = 500

    
    

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


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
    	rundata = setrun()
    prob_data = MultilayerData(rundata)

    rundata.write()
    prob_data.write()
