#!/usr/bin/env python
# encoding: utf-8
""" 
Module to setup multilayer data
""" 

import sys

import numpy as np

from pyclaw import data 

# Multilayer settings
class MultilayerData(data.Data):
    def __init__(self):
        super(MultilayerData,self).__init__()
        
        # Physics parameters
        self.add_attribute('layers',2)
        self.add_attribute('rho',[1025.0,1028.0])
        
        # Algorithm parameters
        self.add_attribute('eigen_method',4)
        self.add_attribute('inundation_method',2)
        self.add_attribute('richardson_tolerance',0.95)
        self.add_attribute('wave_tolerance',[1e-1,2e-1])
        self.add_attribute('dry_limit',False)
    
        # Initial conditions
        self.add_attribute('eta',[0.0,-200.0])
        self.add_attribute('init_type',0)
        self.add_attribute('init_location',[300e3,0.0])
        self.add_attribute('wave_family',4)
        self.add_attribute('epsilon',0.4)
        self.add_attribute('sigma',25e3)
        self.add_attribute('theta',0.0)
        
        # Bathymetry
        self.add_attribute('bathy_type',0)
        self.add_attribute('bathy_location',-50e3)
        self.add_attribute('bathy_left',-4000.0)
        self.add_attribute('bathy_right',-200.0)
        self.add_attribute('bathy_angle',0.0)
        
        self.add_attribute('x0',350e3)
        self.add_attribute('x1',450e3)
        self.add_attribute('x2',480e3)
        self.add_attribute('basin_depth',-3000.0)
        self.add_attribute('shelf_depth',-100.0)
        self.add_attribute('beach_slope',0.008)
        self.add_attribute('h',100)
    
    def write(self,out_file='./multilayer.data',datasource="setrun.py"):
        print "Creating data file %s" % out_file
        
        out_file = open(out_file,'w')
        
        data.data_write(out_file,self,'layers','(Number of layers)')
        data.data_write(out_file,self,'rho','(Densities of layers)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'eigen_method','(Method for calculating eigenspace)')
        data.data_write(out_file,self,'inundation_method','(Method for calculating inundation eigenspace)')
        data.data_write(out_file,self,'richardson_tolerance','(Tolerance for Richardson number)')
        data.data_write(out_file,self,'wave_tolerance','(Tolerance for wave height refinement)')
        data.data_write(out_file,self,'dry_limit','(Turn off limiting when near a dry state)')
        data.data_write(out_file,self,None)
        
        data.data_write(out_file,self,'eta','(Top surface of each layer)')
        data.data_write(out_file,self,'init_type','(Type of initial condition)')
        if self.init_type > 0:
            data.data_write(out_file,self,'epsilon','(Perturbation strength)')
            if self.init_type <= 2 or self.init_type == 5:
                data.data_write(out_file,self,'init_location','(Location for perturbation)')
                data.data_write(out_file,self,'wave_family','(Wave family of the perturbation)')
                if 2 == self.init_type or self.init_type == 5:
                    data.data_write(out_file,self,'angle','(Angle of direction of travel from x-axis)')
                    data.data_write(out_file,self,'sigma','(Gaussian width')
            elif self.init_type == 3:
                data.data_write(out_file,self,'init_location','(Location for perturbation)')
                data.data_write(out_file,self,'sigma','(Gaussian width')
        data.data_write(out_file,self,None)
        
        data.data_write(out_file,self,'bathy_type','(Type of bathymetry prescribed)')
        if self.bathy_type == 0:
            pass
        elif self.bathy_type == 1:
            data.data_write(out_file,self,'bathy_location','(Bathymetry jump location)')
            data.data_write(out_file,self,'bathy_left','(Depth to left of bathy_location)')
            data.data_write(out_file,self,'bathy_right','(Depth to right of bathy_location)')
            data.data_write(out_file,self,'bathy_angle','(Angle from the vertical to put jump)')
        elif self.bathy_type == 2 or self.bathy_type == 3: 
            data.data_write(out_file,self,'x0','(Location of basin end)')
            data.data_write(out_file,self,'x1','(Location of shelf slope end)')
            data.data_write(out_file,self,'x2','(Location of beach slope)')
            data.data_write(out_file,self,'basin_depth','(Depth of basin)')
            data.data_write(out_file,self,'shelf_depth','(Depth of shelf)')
            data.data_write(out_file,self,'beach_slope','(Slope of beach)')
        else:
            print "Invalid bathy_type %s requested, aborting." % self.bathy_type
            out_file.close()
            sys.exit(1)
        
        out_file.close()