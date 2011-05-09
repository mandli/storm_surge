#!/usr/bin/env python
# encoding: utf-8
""" 
Module to setup multilayer data
""" 

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
        
        # Bathymetry
        self.add_attribute('bathy_type',0)
        self.add_attribute('bathy_location',-50e3)
        self.add_attribute('bathy_left',-4000.0)
        self.add_attribute('bathy_right',-200.0)
        
        self.add_attribute('x0',350e3)
        self.add_attribute('x1',450e3)
        self.add_attribute('x2',480e3)
        self.add_attribute('basin_depth',-3000.0)
        self.add_attribute('shelf_depth',-100.0)
        self.add_attribute('beach_slope',0.05)
        self.add_attribute('h',100)
    
    def write(self,out_file='./multilayer.data',datasource="setrun.py"):
        print "Creating data file %s" % out_file
        
        out_file = open(out_file,'w')
        
        data.data_write(out_file,self,'layers','(Number of layers)')
        data.data_write(out_file,self,'rho','(Densities of layers)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'eigen_method','(Method for calculating eigenspace)')
        data.data_write(out_file,self,'richardson_tolerance','(Tolerance for Richardson number)')
        data.data_write(out_file,self,'wave_tolerance','(Tolerance for wave height refinement)')
        data.data_write(out_file,self,'dry_limit','(Turn off limiting when near a dry state)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'eta','(Top surface of each layer)')
        data.data_write(out_file,self,'init_type','(Type of initial condition)')
        data.data_write(out_file,self,'init_location','(Location for perturbation)')
        data.data_write(out_file,self,'wave_family','(Wave family of the perturbation)')
        data.data_write(out_file,self,'epsilon','(Perturbation strength)')
        data.data_write(out_file,self,'sigma','(Gaussian width for init_type=2,3)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'bathy_type','(Type of bathymetry prescribed)')
        if self.bathy_type == 0:
            pass
        elif self.bathy_type == 1:
            data.data_write(out_file,self,'bathy_location','(Bathymetry jump location)')
            data.data_write(out_file,self,'bathy_left','(Depth to left of bathy_location)')
            data.data_write(out_file,self,'bathy_right','(Depth to right of bathy_location)')
        elif self.bathy_type == 2:
            data.data_write(out_file,self,'x0','(Location of basin end)')
            data.data_write(out_file,self,'x1','(Location of shelf slope end)')
            data.data_write(out_file,self,'x2','(Location of beach slope)')
            data.data_write(out_file,self,'basin_depth','(Depth of basin)')
            data.data_write(out_file,self,'shelf_depth','(Depth of shelf)')
            data.data_write(out_file,self,'beach_slope','(Slope of beach)')
            data.data_write(out_file,self,'h','(Height of mid-slope jump near eta(2))')
        else:
            print "Invalid bathy_type %s requested, aborting." % self.bathy_type
            out_file.close()
            sys.exit(1)
        
        out_file.close()