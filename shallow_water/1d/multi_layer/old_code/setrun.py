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

        self.add_attribute('r',0.98)
        self.add_attribute('rho_air',1.15e-3)
        self.add_attribute('eigenvalue_type',2)
        self.add_attribute('dry_tolerance',1e-3)
        
        self.add_attribute('init_type',0)
        self.add_attribute('total_depth',1.0)
        self.add_attribute('depth_ratio',0.75)
        self.add_attribute('init_location',0.5)
        self.add_attribute('beta',0.1)
        self.add_attribute('sigma',0.02)

        self.add_attribute('bathy_type',1)
        self.add_attribute('depth',1.0)

        self.add_attribute('wind_type',2)
        self.add_attribute('A',5.0)
        self.add_attribute('omega',2.0)
        self.add_attribute('N',2.0)
        self.add_attribute('t_length',rundata.clawdata.tfinal-rundata.clawdata.t0)
        self.add_attribute('B',1.5)
        self.add_attribute('Pn',1005.0)
        self.add_attribute('Pc',950.0)
        self.add_attribute('hurricane_velocity',5.0)
        self.add_attribute('R_eye_init',0.0)
        self.add_attribute('ramp_up_time',0.0)

        tf = 10.0
        
    def write(self,file='./problem.data',datasource='setrun.py'):
        """Write out the data file to the path given"""

        print "Creating data file %s" % file
        out_file = data.open_datafile(file)
        
        data.data_write(out_file,self,'r','(Density ratio in layers)')
        data.data_write(out_file,self,'rho_air','(Density of air)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'eigenvalue_type','(Type of eigenvalues to use)')
        data.data_write(out_file,self,'dry_tolerance','(How shallow water is when considered dry)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'init_type','(Type of initial condition)')
        data.data_write(out_file,self,'total_depth','(Total depth of fluid)')
        data.data_write(out_file,self,'depth_ratio','(Ratio of bottom layer to total depth)')
        if self.init_type > 0:
            data.data_write(out_file,self,'init_location','(Location of the initial condition)')
        if self.init_type == 1:
            data.data_write(out_file,self,'beta','(Strength of jump)')
        elif self.init_type == 2:
            data.data_write(out_file,self,'beta','(Height of Gaussian)')
            data.data_write(out_file,self,'sigma','(Width of Gaussian)')
        elif self.init_type == 3:
            data.data_write(out_file,self,'beta','(Strength of perturbation)')
            data.data_write(out_file,self,'wave_family','(Which wave family the perturbation is in)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'bathy_type','(Type of bathymetry to use)')
        if self.bathy_type == 1:
            data.data_write(out_file,self,'depth','(Depth of water, constant)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'wind_type','(Type of wind field to use)')
        if self.wind_type == 1:
            data.data_write(out_file,self,'A','(Wind speed)')
        elif self.wind_type == 2:
            data.data_write(out_file,self,'A','(Hurricane strength parameters)')
            data.data_write(out_file,self,'B','')
            data.data_write(out_file,self,'Pn','(Nominal pressure)')
            data.data_write(out_file,self,'Pc','(Central pressure)')
            self.ramp_up_time = -self.ramp_up_time
            data.data_write(out_file,self,'ramp_up_time','(Time over which the hurricane wind field ramps up to full strength)')
            data.data_write(out_file,self,'hurricane_velocity','(Speed of hurricane)')
            data.data_write(out_file,self,'R_eye_init','(Location of hurricane at t=0)')
        elif self.wind_type == 3:
            data.data_write(out_file,self,'A','(Amplitude of wind)')
            data.data_write(out_file,self,'N','(Number of periods within domain)')
            data.data_write(out_file,self,'omega','(Speed of modulation)')
            data.data_write(out_file,self,'t_length',"(Length of simulation time)")
        
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

    # Number of grid cells:
    clawdata.mx = 100

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

    clawdata.outstyle = 1
    
    # For outstyle 1 specified in hours
    num_hours = 40
    step = 0.25

    if clawdata.outstyle==1:
        # Output nout frames at equally spaced times up to tfinal:
        # clawdata.nout = int(num_hours / step) + int(np.ceil(RAMP_UP_TIME / (step*60**2)))
        # clawdata.tfinal = num_hours * 60.0**2
        clawdata.nout = 160
        clawdata.tfinal = 10.0

    elif clawdata.outstyle == 2:
        # Specify a list of output times.  
        clawdata.tout =  [0.5, 1.0]   # used if outstyle == 2
        clawdata.nout = len(clawdata.tout)

    elif clawdata.outstyle == 3:
        # Output every iout timesteps with a total of ntot time steps:
        iout = 1
        ntot = 60
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
    clawdata.dt_initial = 0.1
    
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
    clawdata.src_split = 1
    
    
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
    
    clawdata.mthbc_xlower = 3
    clawdata.mthbc_xupper = 3
    
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