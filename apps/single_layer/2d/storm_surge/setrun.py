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

import hurricane_data
import multilayer_data
import topo_data

# Ramp up constants
RAMP_UP_TIME = 12*60**2

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
    # Level 1 - 800e3/(80*2) = 5 km
    # Level 2 - 800e3/(80*2*4) = 1.25 km
    # Level 3 - 800e3/(80*2*4*2) = 625 m
    clawdata.xlower = -200e3
    clawdata.xupper = 500e3
    
    clawdata.ylower = -300e3
    clawdata.yupper = 300e3
        

    # Number of grid cells:
    clawdata.mx = 70
    clawdata.my = 60
    # clawdata.mx = 140
    # clawdata.my = 120
    # clawdata.mx = 100
    # clawdata.my = 100

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 8
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.mcapa = 0
    
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = -RAMP_UP_TIME
    
    
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
        clawdata.nout = int(num_hours / step) + int(np.ceil(RAMP_UP_TIME / (step*60**2)))
        clawdata.tfinal = num_hours * 60.0**2
        # # Output nout frames at equally spaced times up to tfinal:
        # clawdata.nout = 48
        # clawdata.tfinal = 60.0**2*clawdata.nout

    elif clawdata.outstyle == 2:
        # Specify a list of output times.
        # t_start = 0.7650e05
        # dt = 0.12E+03 
        # clawdata.tout = [0.0,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]
        # clawdata.tout = [x*60**2 for x in clawdata.tout]
        clawdata.tout = [0.0,1.0,2.0,3.0,4.0]
        # for i in xrange(60):
        #     clawdata.tout.append(i*dt+t_start)
        clawdata.nout = len(clawdata.tout)
    elif clawdata.outstyle == 3:
        # Output every iout timesteps with a total of ntot time steps:
        iout = 1
        ntot = 100
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
    clawdata.dt_initial = 1.e2
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used, and max to allow without 
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.8
    clawdata.cfl_max = 0.9
    
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
    clawdata.mwaves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(mthlim) == mwaves
    clawdata.mthlim = [3,3,3]
    
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
    mxnest = 5

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
    geodata.drytolerance = 1.e-2
    geodata.wavetolerance = 5e-1
    geodata.depthdeep = 2.e2
    geodata.maxleveldeep = 4
    geodata.coeffmanning = 0.025
    # geodata.frictiondepth = 20e1
    geodata.frictiondepth = 1e10
    
    geodata.icoriolis = 1
    geodata.ifriction = 2

    # == settopo.data values ==
    geodata.topofiles = []
    # for topography, append lines of the form
    #   [topotype, minlevel, maxlevel, t1, t2, fname]
    # geodata.topofiles.append([2, 1, 1, 0., 1.e10, 'bowl.topotype2'])
    geodata.topofiles.append([1, 1, 5, 0., 1e10, 'topo.data'])
    
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
    # geodata.regions.append([4, 4, -RAMP_UP_TIME, 1.e10, 400e3,525e3,-325e3,325e3])

    # == setgauges.data values ==
    geodata.gauges = []
    # for gauges append lines of the form  [gaugeno, x, y, tstart, tend]
    N_gauges = 21
    for i in xrange(0,N_gauges):
        x = 455e3 # This is right where the shelf turns into beach, 100 meter of water
        # x = -80.0 * (23e3 / 180) + 500e3 - 5e3  # 1 km off shore
        y = 550e3 / (N_gauges + 1) * (i+1) + -275e3       # Start 25 km inside domain
        geodata.gauges.append([i, x, y, 0.0, 1e10])
        print "Gauge %s: (%s,%s)" % (i,x/1e3,y/1e3)

    # == setfixedgrids.data values ==
    geodata.fixedgrids = []
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]
    # geodata.fixedgrids.append([1., 2., 4, 0., 100., 0., 100., 11, 11, 0, 0])
    
    return rundata
    # end of function setgeo
    # ----------------------
        
def set_hurricane_data(ramp_up_time=RAMP_UP_TIME):
    data = hurricane_data.HurricaneData(ramp_up_time)
    
    # Source terms to be included
    data.wind_src = True
    data.pressure_src = True
    
    # Momentum based refinement
    data.momentum_refinement = False
    data.max_speed_nest = 5
    # data.speed_nest = [1.0,2.0,3.0,4.0,5.0]
    # data.speed_nest = [0.5,1.0,2.0,3.0,4.0,5.0]
    data.speed_nest = [0.25,0.5,1.0,2.0,3.0,4.0]
    
    # Hurricane location based refinement
    data.max_R_nest = 4
    data.R_refine = [60.0e3,50e3,40e3,30e3]
        
    # Wind strength based refinement
    data.max_wind_nest = 3
    # data.wind_refine = [0.001,0.005,0.001]
    data.wind_refine = [20.0,40.0,60.0]
    
    # Pressure source term tolerance for gradient value
    data.pressure_tolerance = 1e-4
    
    # Ramp up time for hurricane
    data.ramp_up_t = RAMP_UP_TIME
    
    # Path of hurricane, speed in m/s
    velocity = 5.0
    angle = 0.0
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
    data.rho_air = 1.15     # Density of air, this also includes 
    
    return data

def set_multilayer_data():
    data = multilayer_data.MultilayerData()
    
    # Physical parameters
    data.layers = 1
    data.rho = [1025.0]
    # The rest of these are ignored for single layers
    
    # Algorithm Parameters
    data.eigen_method = 1
    data.richardson_tolerance = 0.95
    data.wave_tolerance = 1e-1
    
    # Initial conditions
    data.eta = [0.0]
    data.init_type = 0
    data.init_location = [300e3,0.0]
    data.wave_family = 4
    data.epsilon = 0.4
    data.sigma = 25e3
    
    # Bathy settings
    data.bathy_type = 0
    
    return data

# =====================
#  Main run function 
# =====================
if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
        rundata = setrun()
    hurricane_data = set_hurricane_data()
    multilayer_data = set_multilayer_data()

    # Write our run data file
    rundata.write()
    hurricane_data.write()
    multilayer_data.write()
    
    # Write out topography and qinit data files if needed
    topo_file = './topo.data'
    topo_data.write_topo_file(topo_file,topo_type=1,factor=4,
                            bathy_type='simple_shelf',plot=True,force=True)
    
