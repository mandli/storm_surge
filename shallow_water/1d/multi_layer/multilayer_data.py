#!/usr/bin/env python
# encoding: utf-8
r"""
Data class for multilayer runs in 1d

:Authors:
	Kyle Mandli (2011-03-18) Initial version
"""
# ============================================================================
#      Copyright (C) 2011 Kyle Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import sys
import os

import pyclaw.data as data

# Simple hurricane data format
class MultilayerData(data.Data):
    def __init__(self):
        super(MultilayerData,self).__init__()
        
        # Physics
        self.add_attribute('rho_air',1.0)
        self.add_attribute('rho_1',1025.0)
        # self.add_attribute('rho_2',self.rho_1/0.99)
        self.add_attribute('rho_2',1028.0)
        
        # Algorithm
        self.add_attribute('dry_tolerance',1e-3)
        self.add_attribute('eigen_method',4)
        self.add_attribute('inundation_method',1)
        
        # Initial condition
        self.add_attribute('init_type',4)
        self.add_attribute('init_location',300e3)
        self.add_attribute('wave_family',4)
        self.add_attribute('eta_1',0.0)
        self.add_attribute('eta_2',-300)
        # self.add_attribute('epsilon',0.027)
        self.add_attribute('epsilon',0.4)
        self.add_attribute('sigma',25e3)
        
        # Bathymetry
        self.add_attribute('bathy_location',-30e3)
        self.add_attribute('bathy_left',-4000.0)
        self.add_attribute('bathy_right',-200.0)
        
        # Wind
        self.add_attribute('wind_type',0)
        self.add_attribute('A',5.0)
        self.add_attribute('omega',2.0)
        self.add_attribute('N',2.0)
        self.add_attribute('t_length',10.0)
        self.add_attribute('B',1.5)
        self.add_attribute('Pn',1005.0)
        self.add_attribute('Pc',950.0)
        self.add_attribute('hurricane_velocity',5.0)
        self.add_attribute('R_eye_init',0.0)
        self.add_attribute('ramp_up_time',0.0)
        
    def write(self,file='./problem.data',datasource='setrun.py'):
        """Write out the data file to the path given"""

        print "Creating data file %s" % file
        out_file = data.open_datafile(file)
        
        data.data_write(out_file,self,'rho_air','(Density of air)')
        data.data_write(out_file,self,'rho_1','(Density of top layer)')
        data.data_write(out_file,self,'rho_2','(Density of bottom layer)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'dry_tolerance','(Dry state tolerance)')
        data.data_write(out_file,self,'eigen_method','(Method for calculating eigenspace)')
        data.data_write(out_file,self,'inundation_method','(Method for calculating inundation Riemann problem)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'init_type','(Type of initial condition)')
        data.data_write(out_file,self,'init_location','(Location for perturbation)')
        data.data_write(out_file,self,'wave_family','(Wave family of the perturbation)')
        data.data_write(out_file,self,'eta_1','(Steady state top surface)')
        data.data_write(out_file,self,'eta_2','(Steady state internal surface)')
        data.data_write(out_file,self,'epsilon','(Perturbation strength)')
        data.data_write(out_file,self,'sigma','(Gaussian width for init_type=2,3)')
        data.data_write(out_file,self,None)
        data.data_write(out_file,self,'bathy_location','(Bathymetry jump location)')
        data.data_write(out_file,self,'bathy_left','(Depth to left of bathy_location)')
        data.data_write(out_file,self,'bathy_right','(Depth to right of bathy_location)')
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
