#!/usr/bin/env python
# encoding: utf-8
""" 
Functions for creating topography
    
"""

import os
import numpy as np

from pyclaw.data import Data
from pyclaw.util import create_topo_func
from pyclaw import util
import pyclaw.geotools.topotools as tt

# ============================================================================
#  Topography generation functions 
# ============================================================================
def write_topo_file(topo_file,bathy_type=None,plot=False,force=False,verbose=False):
    """Creates topography file needed by the simulation"""
    
    if os.path.exists(topo_file):
        if force:
            ans = raw_input("The file %s already exists, replace: [Y/n] " % topo_file)
            if ans[0].capitalize() == 'N':
                print "Exiting"
                return
        else:
            print "File %s already exists, use force if you want to replace the existing file." % topo_file
    print "Creating topography file ",topo_file
    
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
    
    # ========================================================================
    # Bathy parameters
    beach_slope = 0.008
    basin_depth = -3000
    shelf_depth = -200
    x0 = 350e3
    x1 = 450e3
    x2 = 480e3
    
    shelf_slope = (basin_depth - shelf_depth) / (x0 - x1)
    y_end = beach_slope * (xupper - 477e3) - 100.0
    
    # Points and depths for gulf shelf
    #  1: (25°39'2.85"N, 86° 7'24.77"W)   --   -3228 m   --   0.0 m
    #  2: (27°53'44.74"N, 88° 0'34.02"W)   --   -2438 m   --   312.17313 km
    #  3: (28°59'47.14"N, 88°59'53.19"W)   --   -188 m   --    467.59957 km
    #  4: ( 29° 4'6.90"N,  89° 4'11.39"W)    --   0 m   --   479.10557 km
    
    bathy_profiles = {"simple_shelf":[(xlower,basin_depth),
                                      (x0,basin_depth),
                                      (x1,shelf_depth),
                                      (x2,shelf_depth),
                                      (xupper,beach_slope*(xupper-x2)+shelf_depth)],
                      "shallow_shelf":[(477e3,-100),
                                       (xupper,y_end)],
                      "gulf_shelf":[(0.0,-3228),
                                    (312e3,-2438),
                                    (467e3,-188),
                                    (479e3,0.0),
                                    (579e3,300.0)],
                      "step_shelf1":[(0.0,-2000.0),
                                     (470e3-0.001,-2000.0),
                                     (470e3,-200.0),
                                     (500e3,-200.0)],
                      "shelf":[(0.0,-3000),
                               (400e3,-2700),
                               (450e3,-100),
                               (500e3,-100)],
                      "continental_shelf":[(2000e3,-7000),
                                           (2800e3,-3000),
                                           (2900e3,-100),
                                           (3000e3,0.0)],
                      "flat":[(0.0,-2000),
                              (400e3,-2000)]}

    # Pick out bathymetry profile
    if bathy_type is None:
        bathy_profile = 'flat'
    if bathy_type in bathy_profiles.keys():
        bathy_profile = bathy_profiles[bathy_type]
    if verbose:
        print "%s profile: %s" % (bathy_type,bathy_profile)

    # Create and write topography
    bathy_func = util.create_topo_func(bathy_profile,verbose=verbose)
    N = len(bathy_profile)
    factor = 4
    tt.topo1writer(topo_file,bathy_func,xlower,xupper,ylower,yupper,mx,my)

    # Plotting
    if plot:
        import matplotlib.pyplot as plt
        (X,Y,Z) = tt.topofile2griddata(topo_file,topotype=topo_type)
        slice_index = int(X.shape[0] / 2.0)
        x = X[slice_index,:]
        b = Z[slice_index,:]
        
        # Plotting parameters
        limits = [np.min(x),np.max(x),np.min(b)*1.05,np.max(b)*1.05]
        wet_index = np.nonzero(b < 0.0)
        
        # Profile plot
        plt.figure(1)
        plt.hold(True)
        plt.plot(x,b,'k',linewidth=2)
        plt.plot(x[wet_index],np.zeros(x[wet_index].shape),'b')
        plt.axis(limits)
        # 
        # plt.plot(x,bath[index,:].T,'k',linewidth=2)
        # for (i,surface) in enumerate(eta):
        #     surface_vector = np.ones(bath.shape[1]) * surface
        #     wet_index = np.nonzero(bath[index,:] <= surface_vector)
        #     plt.plot(x[wet_index],surface_vector[wet_index],'k',linewidth=1)            
        # plt.fill_between(x,bath[index,:].T,x*eta[0],where=(bath[index,:].T < 0.0),color=[0,0.5,1])
        # plt.axis([-200e3,500e3,ylimits[0],ylimits[1]])
        plt.title("Bathymetry Cross-Section")
        plt.xlabel('km')
        plt.ylabel('m')
        locs,labels = plt.xticks()
        labels = locs/1e3
        plt.xticks(locs,labels)
        plt.hold(False)
        plt.show()
        
if __name__ == "__main__":
    write_bathy_profile('./topo.data',None,plot=True,force=True)