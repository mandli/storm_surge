#!/usr/bin/env python
# encoding: utf-8
""" 
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
""" 

import os
import numpy as np

from pyclaw import data 
import pyclaw.util as util

import multilayer_data
import hurricane_data
import topo_data

direction = 2

#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------
    
    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """ 
    
    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    ndim = 2
    rundata = data.ClawRunData(claw_pkg, ndim)

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------

    rundata = setgeo(rundata)   # Defined below     
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.ndim = ndim
    
    # Lower and upper edge of computational domain:
    clawdata.xlower = 0.0
    clawdata.xupper = 1.0
    
    clawdata.ylower = 0.0
    clawdata.yupper = 1.0
        

    # Number of grid cells:
    # clawdata.mx = 70
    factor = 1
    if direction == 1:
        clawdata.mx = 70 * factor
        clawdata.my = 60
    elif direction == 2:
        clawdata.mx = 70
        clawdata.my = 60 * factor
    # clawdata.my = 120
    # clawdata.mx = 100
    # clawdata.my = 100

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 6

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 8
    
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
    # Number of hours to simulate
    num_hours = 40
    # Output interval per hour, 1 = every hour, 0.5 = every half hour, etc...
    step = 0.25
    
    if clawdata.outstyle==1:
        clawdata.nout = 80
        # if wave_family < 4 and wave_family > 1:
        clawdata.tfinal = 1.0
        # else:
            # clawdata.tfinal = 0.1

    elif clawdata.outstyle == 2:
        # Specify a list of output times.
        t_start = 0.71100e5
        dt = 0.4046e2
        # clawdata.tout = [0.0,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]
        # clawdata.tout = [x*60**2 for x in clawdata.tout]
        clawdata.tout = []
        for i in xrange(100):
            clawdata.tout.append(i*dt+t_start)
        clawdata.nout = len(clawdata.tout)
    elif clawdata.outstyle == 3:
        # Output every iout timesteps with a total of ntot time steps:
        iout = 1
        ntot = 80
        clawdata.iout = [iout, ntot]
    


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 2
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = 1
    
    # Initial time step for variable dt.  
    # If dt_variable==0 then dt=dt_initial for all steps:
    # clawdata.dt_initial = 0.64e2
    clawdata.dt_initial = 0.575e-03
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used, and max to allow without 
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.8
    clawdata.cfl_max = 0.9
    # clawdata.cfl_desired = 0.4
    # clawdata.cfl_max = 0.5
    
    # Maximum number of time steps to allow between output times:
    clawdata.max_steps = 5000

    
    

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Transverse order for 2d or 3d (not used in 1d):
    clawdata.order_trans = 2
    
    # Number of waves in the Riemann solution:
    clawdata.mwaves = 6
    
    # List of limiters to use for each wave family:  
    # Required:  len(mthlim) == mwaves
    clawdata.mthlim = [3,3,3,3,3,3]
    
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
    
    clawdata.mthbc_xlower = 1
    clawdata.mthbc_xupper = 1
    
    clawdata.mthbc_ylower = 1
    clawdata.mthbc_yupper = 1
    

    # ---------------
    # AMR parameters:
    # ---------------


    # max number of refinement levels:
    mxnest = 1

    clawdata.mxnest = -mxnest   # negative ==> anisotropic refinement in x,y,t

    # List of refinement ratios at each level (length at least mxnest+1)
    clawdata.inratx = [2,2,2,2,2]
    clawdata.inraty = [2,2,2,2,2]
    clawdata.inratt = [2,2,2,2,2]
    # clawdata.inratt = [1,1,1]


    # Specify type of each aux variable in clawdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    clawdata.auxtype = ['center','center','center','center','center','center',
                        'center','center']


    clawdata.tol = -1.0     # negative ==> don't use Richardson estimator
    clawdata.tolsp = 0.5    # used in default flag2refine subroutine
                            # (Not used in geoclaw!)

    clawdata.cutoff = 0.7   # efficiency cutoff for grid generation
    clawdata.kcheck = 2     # how often to regrid (every kcheck steps)
    clawdata.ibuff  = 2     # width of buffer zone around flagged points

    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    return rundata
    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geodata = rundata.geodata
    except:
        print "*** Error, this rundata has no geodata attribute"
        raise AttributeError("Missing geodata attribute")

    # == setgeo.data values ==
    geodata.igravity = 1
    geodata.gravity = 9.81
    geodata.icoordsys = 1

    # == settsunami.data values ==
    geodata.sealevel = 0.
    # geodata.drytolerance = 1.e-3
    geodata.drytolerance = 1.e-3
    geodata.wavetolerance = 5e-1
    # geodata.wavetolerance = 1e0
    geodata.depthdeep = 2.e2
    geodata.maxleveldeep = 4
    geodata.ifriction = 1
    geodata.coeffmanning = 0.025
    geodata.frictiondepth = 20.0

    # == settopo.data values ==
    geodata.topofiles = []
    # for topography, append lines of the form
    #   [topotype, minlevel, maxlevel, t1, t2, fname]
    # geodata.topofiles.append([2, 1, 1, 0., 1.e10, 'bowl.topotype2'])
    # geodata.topofiles.append([1, 1, 5, 0., 1e10, 'topo.data'])
    
    # == setdtopo.data values ==
    geodata.dtopofiles = []
    # for moving topography, append lines of the form:  (<= 1 allowed for now!)
    #   [minlevel,maxlevel,fname]

    # == setqinit.data values ==
    geodata.qinitfiles = []  
    geodata.iqinit = 0
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]
    # geodata.qinitfiles.append([1,3,'qinit.data'])

    # == setregions.data values ==
    geodata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    # geodata.regions.append([1, 1, 0.e0, 1.e10, -100.,100., -100.,100.])
    # geodata.regions.append([1, 2, 0.e0, 1.e10,    0.,100., -100.,100.])

    # == setgauges.data values ==
    geodata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, tstart, tend]
    # N_gauges = 21
    # for i in xrange(0,N_gauges):
    #     x = -80.0 * (23e3 / 180) + 500e3 - 5e3  # 1 km off shore
    #     y = 550e3 / (N_gauges + 1) * (i+1) + -275e3       # Start 25 km inside domain
    #     geodata.gauges.append([i, x, y, 0.0, 1e10])
    #     print "Gauge %s: (%s,%s)" % (i,x/1e3,y/1e3)

    # == setfixedgrids.data values ==
    geodata.fixedgrids = []
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]
    # geodata.fixedgrids.append([1., 2., 4, 0., 100., 0., 100., 11, 11, 0, 0])
    
    return rundata
    # end of function setgeo
    # ----------------------
        
def set_hurricane_data(ramp_up_time=0.0):
    
    data = hurricane_data.HurricaneData(0.0)
    
    # Source terms to be included
    data.wind_src = False
    data.pressure_src = False
    
    # Momentum based refinement
    data.momentum_refinement = False
    data.max_speed_nest = 5
    # data.speed_nest = [1.0,2.0,3.0,4.0,5.0]
    # data.speed_nest = [0.5,1.0,2.0,3.0,4.0,5.0]
    data.speed_nest = [0.25,0.5,1.0,2.0,3.0,4.0]
    
    # Hurricane location based refinement
    data.max_R_nest = 3
    data.R_refine = [60.0e3,40e3,20e3]
        
    # Wind strength based refinement
    data.max_wind_nest = 0
    # data.wind_refine = [0.001,0.005,0.001]
    data.wind_refine = [20.0,40.0,60.0]
    
    # Pressure source term tolerance
    data.pressure_tolerance = 1e-4
    
    # Ramp up time for hurricane
    data.ramp_up_t = 0.0
    
    # Type of wind field
    data.wind_type = 1
    data.wind_tolerance = 1e-6
    
    # Path of hurricane, speed in m/s
    velocity = 5.0
    angle = 0.0 * np.pi
    # Speeds of hurricane
    data.hurricane_velocity = (velocity * np.cos(angle),velocity * np.sin(angle)) 
    # Initial position of hurricane eye at t = 0
    data.R_eye_init = (0.0,0.0) 

    # Hurricane parameters
    # These match Hurricane Tracy
    data.A = 23.0           # Hurricane model fit parameter
    data.B = 1.5     
    data.Pn = 1005.0        # Nominal atmospheric pressure     
    data.Pc = 950.0         # Pressure in the eye of the hurricane    
    data.rho_air = 1.15     # Density of air
    
    return data
    
def set_multilayer_data():
    data = multilayer_data.MultilayerData()
    
    # Physical parameters
    data.layers = 2
    # data.rho = [1.0,0.0]
    data.rho = [0.90,1.0]
    
    # Algorithm Parameters
    data.eigen_method = 4
    data.richardson_tolerance = 0.95
    data.wave_tolerance = [0.1,0.1]
    data.dry_limit = True
    
    # Initial conditions
    # data.eta = [0.0,-5000.0]
    data.eta = [0.0,-0.6]
    if direction == 1:
        data.init_type = 2
    elif direction == 2:
        data.init_type = 6
    data.init_location = [0.25,0.25]
    data.wave_family = 0
    data.epsilon = 0.02
    data.sigma = 0.02
    
    # Bathy settings
    data.bathy_location = 0.5
    data.bathy_left = -1.0
    data.bathy_right = -0.2
    
    return data
    
# =====================
#  Main run function 
# =====================
if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    # if len(sys.argv) == 2:
    #     rundata = setrun(sys.argv[1])
    # else:
    rundata = setrun()
    hurricane_data = set_hurricane_data()
    multilayer_data = set_multilayer_data()

    # Write our run data file
    rundata.write()
    hurricane_data.write()
    multilayer_data.write()    
    
    # Write out topography and qinit data files if needed
    topo_data.write_topo_file('./topo.data',bathy_type='flat',
                                        plot=False,force=False)
    
