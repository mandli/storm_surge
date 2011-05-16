#!/usr/bin/env python
# encoding: utf-8
""" 
Functions for creating topography
    
"""

import os

import numpy as np

from pyclaw import data
from pyclaw import util

# ============================================================================
#  Topography generation functions 
# ============================================================================
def write_topo_file(file,bathy_type=None,plot=False,force=False):
    """Creates topography file needed by the simulation"""

    from pyclaw.data import Data
    from pyclaw.util import create_topo_func

    if os.path.exists(file) and force:
        ans = raw_input("The file %s already exists, replace: [Y/n] ")
        if ans[0].capitalize() == 'N':
            print "Exiting"
            return
    print "Creating topography file ",file
    
    # Parameters
    data = Data('./amr2ez.data')
    dx = abs(data.xupper-data.xlower) / (data.mx)
    dy = abs(data.yupper-data.ylower) / (data.my)
    d = min(dx,dy)
    mx = int((data.xupper-data.xlower) / d) + 8
    my = int((data.yupper-data.ylower) / d) + 8
    
    xlower = data.xlower-d*4.0
    ylower = data.ylower-d*4.0
    xupper = data.xupper+d*4.0
    yupper = data.yupper+d*4.0
    
    topo_type = 1
    
    # Bathy types
    # Beach slopes
    beach_slope = 0.05
    y_end = beach_slope * (xupper - 477e3) - 100.0
    
    # New bathy support
    epsilon = 1.0
    x0 = 350e3
    x1 = 450e3
    x2 = 480e3
    x0 = 350e3
    x1 = 450e3
    x2 = 480e3
    basin_depth = -3000
    shelf_depth = -100
    eta = [0.0,-300]
    h = 100.0
    A = basin_depth - eta[1] + 0.5*h
    B = shelf_depth - eta[1] - 0.5*h
    eta_int = (A*x1 - B*x0) / (A-B)
    shelf_slope =  A / (x0 - eta_int)
    
    # Points and depths for gulf shelf
    #  1: (25°39'2.85"N, 86° 7'24.77"W)   --   -3228 m   --   0.0 m
    #  2: (27°53'44.74"N, 88° 0'34.02"W)   --   -2438 m   --   312.17313 km
    #  3: (28°59'47.14"N, 88°59'53.19"W)   --   -188 m   --    467.59957 km
    #  4: ( 29° 4'6.90"N,  89° 4'11.39"W)    --   0 m   --   479.10557 km
    
    
    bathys = {"new_bathy1":[(-250e3,basin_depth),(x0,basin_depth),
                            (x0,basin_depth),
                            (x1,shelf_depth),
                            (x2,shelf_depth),
                            (500e3,beach_slope*(500e3-x2)+shelf_depth)],
              "new_bathy2":[(-250e3,basin_depth),(x0,basin_depth),
                       (eta_int-epsilon,shelf_slope*(eta_int-epsilon-x0)+basin_depth),
                       (eta_int+epsilon,shelf_depth),
                       (x2,shelf_depth),
                       (500e3,beach_slope*(500e3-x2)+shelf_depth)],
              "new_bathy3":[(-250e3,basin_depth),(x0,basin_depth),
                       (eta_int-epsilon,shelf_slope*(eta_int-epsilon-x0)+basin_depth),
                       (eta_int+epsilon,shelf_slope*(eta_int+epsilon-x1)+shelf_depth),
                       (x1,shelf_depth),(x2,shelf_depth),(500e3,beach_slope*(500e3-x2)+shelf_depth)],
              "shallow_shelf":[(477e3,-100),(xupper,y_end)],
              "gulf_shelf":[(0.0,-3228),(312e3,-2438),(467e3,-188),(479e3,0.0),(579e3,300.0)],
              "step_shelf1":[(0.0,-2000.0),(470e3-0.001,-2000.0),(470e3,-200.0),(500e3,-200.0)],
              "shelf":[(0.0,-3000),(400e3,-2700),(450e3,-100),(500e3,-100)],
              "step_shelf2":[(0.0,-4000.0),(450e3-0.001,-4000.0),(450e3,-200.0),(500e3,-200.0)],
              "continental_shelf":[(2000e3,-7000),(2800e3,-3000),(2900e3,-100),(3000e3,0.0)],
              "flat":[(0.0,-2000),(400e3,-2000)]}
    
    if bathy_type is None:
        bathy_profile = flat
    else:
        if bathy_type in bathys.keys():
            bathy_profile = bathys[bathy_type]
    
    # Create and write topography
    import pyclaw.geotools.topotools as tt
    bathy_func = util.create_topo_func(bathy_profile)
    N = len(bathy_profile)
    tt.topo1writer(file,bathy_func,xlower,xupper,ylower,yupper,N,N)

    # Plotting
    if plot:
        import matplotlib.pyplot as plt
        x = np.linspace(xlower,xupper,mx)
        y = np.linspace(ylower,yupper,my)
        [X,Y] = np.meshgrid(x,y)
        
        # Determine limits
        ylimits = [0.0,0.0]
        for i in xrange(len(bathy_profile)):
            coord = bathy_profile[i][1]
            ylimits = [min(coord,ylimits[0]),max(coord,ylimits[1])]
        ylimits = [ylimits[0]*0.1+ylimits[0],ylimits[1]*0.1+ylimits[1]]
        if ylimits[1] == 0.0:
            ylimits[1] = ylimits[1] + abs(ylimits[0]*0.1)

        bath = np.empty(X.shape)
        bath = bathy_func(X,Y)
        
        # Profile plot
        plt.figure(1)
        index = np.trunc(my/2.0)
        plt.hold(True)
        plt.plot(x,bath[index,:].T,'k',linewidth=2)
        plt.fill_between(x,bath[index,:].T,x*0.0,where=(bath[index,:].T < 0.0),color=[0,0.5,1])
        plt.axis([-200e3,500e3,ylimits[0],ylimits[1]])
        plt.title("Bathymetry Cross-Section")
        plt.xlabel('km')
        plt.ylabel('m')
        locs,labels = plt.xticks()
        labels = locs/1e3
        plt.xticks(locs,labels)
        plt.hold(False)
        plt.show()
        