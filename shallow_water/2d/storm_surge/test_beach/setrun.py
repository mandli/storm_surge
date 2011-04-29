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

# Ramp up constants
RAMP_UP_TIME = 12*60**2

# Simple hurricane data format
class HurricaneData(data.Data):
    def __init__(self):
        super(HurricaneData,self).__init__()

        # Source terms to be included
        self.add_attribute('wind_src',True)
        self.add_attribute('pressure_src',False)
        
        # Momentum based refinement
        self.add_attribute('momentum_refinement',False)
        self.add_attribute('max_speed_nest',5)
        # self.add_attribute('speed_nest',[1.0,2.0,3.0,4.0,5.0])
        self.add_attribute('speed_nest',[0.5,1.0,2.0,3.0,4.0,5.0]) # Just added another level of refinement
        
        # Hurricane location based refinement
        self.add_attribute('max_R_nest',3)
        self.add_attribute('R_refine',[60.0e3,40e3,20e3])
        
        # Wind strength based refinement
        self.add_attribute('max_wind_nest',0)
        # self.add_attribute('wind_refine',[0.001,0.005,0.001])
        self.add_attribute('wind_refine',[20.0,40.0,60.0])
        
        # Pressure source term tolerance
        self.add_attribute('pressure_tolerance',1e-4)
        
        # Ramp up time for hurricane
        self.add_attribute("ramp_up_t",RAMP_UP_TIME)
        
        # Path of hurricane, speed in m/s
        self.add_attribute('hurricane_velocity_X',5.0)  # Speed of hurricane
        self.add_attribute('hurricane_velocity_Y',0.0)
        
        # Initial position of hurricane eye at t = 0
        self.add_attribute('R_eye_init_X',0.0)     # Initial position
        self.add_attribute('R_eye_init_Y',0.0)

        # Hurricane parameters
        # These match Hurricane Tracy
        self.add_attribute('A',23.0)        # Hurricane model fit parameter
        self.add_attribute('B',1.5)     
        self.add_attribute('Pn',1005.0)     # Nominal atmospheric pressure     
        self.add_attribute('Pc',950.0)      # Pressure in the eye of the hurricane    
        self.add_attribute('rho_air',1.15e-3) # Density of air
        
    def write(self,file='./hurricane.data'):
        """Write out the data file to the path given"""

        print "Creating data file %s" % file
        super(HurricaneData,self).write(supplementary_file=file)

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
    # clawdata.mx = 100
    # clawdata.my = 100

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = 9
    
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
        iout = 100
        ntot = 10000
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
    mxnest = 4

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
                        'center','center','center']


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
    geodata.ifriction = 1
    geodata.coeffmanning = 0.025
    geodata.frictiondepth = 20.0

    # == settopo.data values ==
    geodata.topofiles = []
    # for topography, append lines of the form
    #   [topotype, minlevel, maxlevel, t1, t2, fname]
    # geodata.topofiles.append([2, 1, 1, 0., 1.e10, 'bowl.topotype2'])
    geodata.topofiles.append([1, 1, 3, 0., 1e10, 'topo.data'])
    
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
    N_gauges = 21
    for i in xrange(0,N_gauges):
        x = -80.0 * (23e3 / 180) + 500e3 - 5e3  # 1 km off shore
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

# ============================================================================
#  Topography generation functions 
# ============================================================================
def write_topo_file(file,plot=False):
    """Creates topography file needed by the simulation"""

    from pyclaw.data import Data
    from pyclaw.util import create_topo_func

    print "Creating topography file ",file
    
    # Parameters
    data = Data('amr2ez.data')
    dx = abs(data.xupper-data.xlower) / (data.mx * 2**4)
    dy = abs(data.yupper-data.ylower) / (data.my * 2**4)
    d = min(dx,dy)
    mx = int((data.xupper-data.xlower) / d) + 4
    my = int((data.yupper-data.ylower) / d) + 4
    
    xlower = data.xlower-d*2.0
    ylower = data.ylower-d*2.0
    xupper = data.xupper+d*2.0
    yupper = data.yupper+d*2.0
    
    topo_type = 1
    
    # Bathy types
    beach_slope = 0.05
    y_end = beach_slope * (xupper - 477e3) - 100.0
    shallow_shelf = [(477e3,-100),(xupper,y_end)]
    # Points and depths
    #  1: (25°39'2.85"N, 86° 7'24.77"W)   --   -3228 m   --   0.0 m
    #  2: (27°53'44.74"N, 88° 0'34.02"W)   --   -2438 m   --   312.17313 km
    #  3: (28°59'47.14"N, 88°59'53.19"W)   --   -188 m   --    467.59957 km
    #  4: ( 29° 4'6.90"N,  89° 4'11.39"W)    --   0 m   --   479.10557 km
    gulf_shelf = [(0.0,-3228),(312e3,-2438),(467e3,-188),(479e3,0.0),(579e3,300.0)]
    continental_shelf_2 = [(2000e3,-7000),(2800e3,-3000),(2900e3,-100),(3000e3,0.0)]
    flat = [(0.0,-100)]
    
    import pyclaw.geotools.topotools as tt
    bathy_func = util.create_topo_func(shallow_shelf)
    N = len(shallow_shelf)
    tt.topo1writer(file,bathy_func,xlower,xupper,ylower,yupper,N,N)
    # 
    # topo_type = 1
    # out_file = open(file,'w')
    # for point in flat:
    #     out_file.write("%f %f %f\n" % (point[0],upper_coord[1],point[1]))
    #     out_file.write("%f %f %f\n" % (point[0],lower_coord[1],point[1]))
    # out_file.close()
    
    
    # util.write_topo_file(file,bathy_func,mx,my,lower_coord,upper_coord,topo_type)
    
    if plot:
        import matplotlib.pyplot as plot
        # plot.subplot(2,1,1)
        # [x,y,Z] = tt.topofile2griddata(file,topo_type=topo_type)
        # [X,Y] = np.meshgrid(x,y)
        # plot.pcolor(X,Y,Z.T)
        # plot.colorbar()
        # plot.axis('equal')
        # plot.xlabel('km')
        # plot.ylabel('km')
        # locs,labels = plot.xticks()
        # labels = locs/1.e3
        # plot.xticks(locs,labels)
        # locs,labels = plot.yticks()
        # labels = locs/1.e3
        # plot.yticks(locs,labels)
        # plot.subplot(2,1,2)
        plot.hold(True)
        plot.plot(x,Z.T[10,:],'k',x,Z.T[10,:],'ro')
        plot.fill_between(x,Z.T[10,:],where=Z.T[10,:]<0.0,color='b')
        # plot.plot(x,np.zeros(x.shape),'b-')
        plot.axis([data.xlower,data.xupper,-120,100])
        plot.title('Topography Cross-Section')
        plot.xlabel('km')
        plot.ylabel('m')
        locs,labels = plot.xticks()
        labels = locs/1.e3
        plot.xticks(locs,labels)
        plot.show()
        

# =========================
#  Qinit data generation 
# =========================
def gaussian(x,y): 
    """Simple gaussian hump""" 
    center = (0.0,0.0) 
    sigma = (20e3,20e3) 
    # sigma = (0.5,0.5)
    amplitude = 1.0 
     
     
    return amplitude * np.exp( -((x-center[0])/sigma[0])**2 ) \
                     * np.exp( -((y-center[1])/sigma[1])**2 ) 
 
 
def write_qinit_file(file,plot=False): 
    """Creates each of the files needed by the simulation.""" 
    
    from pyclaw.data import Data
    from pyclaw.util import create_topo_func
     
    print "Creating qinit file ", file 
    
    # Parameters
    data = Data('amr2ez.data')
    dx = abs(data.xupper-data.xlower) / float(data.mx)
    dy = abs(data.yupper-data.ylower) / float(data.my)
    d = min(dx,dy)
    mx = int((data.xupper-data.xlower) / d) + 20
    my = int((data.yupper-data.ylower) / d) + 20
    
    lower_coord = (data.xlower-d*10.0,data.ylower-d*10.0)
    upper_coord = (data.xupper+d*10.0,data.yupper+d*10.0)
    
    # Calculate function and write to file 
    eta = np.zeros((mx,my)) 
    qinit_file = open(file,'w') 
    for (j,y) in enumerate(-np.linspace(-upper_coord[1],-lower_coord[1],my)): 
        for (i,x) in enumerate(np.linspace(lower_coord[0],upper_coord[0],mx)): 
            eta[i,j] = gaussian(x,y) 
            qinit_file.write("%s %s %s\n" % (x,y,eta[i,j])) 
    qinit_file.close() 
     
    # Plot qinit if requested 
    if plot: 
        import matplotlib.pyplot as plot 
         
        [X,Y] = np.meshgrid(np.linspace(lower_coord[0],upper_coord[0],mx), 
                            np.linspace(lower_coord[1],upper_coord[1],my)) 
        plot.pcolor(X,Y,eta.T) 
        plot.title("Initial Condition Data") 
        plot.axis("equal") 
        plot.colorbar() 
        plot.show() 

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

    # Write our run data file
    rundata.write()
    
    # Write out hurricane data info
    hurricane_data = HurricaneData()
    hurricane_data.write()
    
    # Write out topography and qinit data files if needed
    topo_file = './topo.data'
    qinit_file = './qinit.data'
    if not(os.path.exists(topo_file)):
        write_topo_file(topo_file,plot=True)
    # if not(os.path.exists(qinit_file)):
    #     write_qinit_file(qinit_file,plot=True)
    
