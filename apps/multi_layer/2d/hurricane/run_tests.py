#!/usr/bin/env python
# encoding: utf-8
r"""
Runs tests for angles relative to a continental shelf

:Authors:
    Kyle Mandli (2011-04-20) Initial version
"""
# ============================================================================
#      Copyright (C) 2011 Kyle Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import os
import sys

import numpy as np

import test_runs

class TwoLayerBaseTest(test_runs.TestML2D):
    
    def __init__(self,rho_ratio=0.97,velocity=5.0,angle=0.0,eye=(0.0,0.0),mxnest=1,
                    mx=70,my=60,bathy_type='simple_shelf',friction=True,eta=[0.0,-300.0]):
                    
        super(TwoLayerBaseTest,self).__init__()
        
        self.type = "storm_surge"
        self.name = "multi_layer"
        self.setplot = "setplot"
        
        self.run_data.clawdata.mx = mx
        self.run_data.clawdata.my = my
        self.run_data.clawdata.mxnest = -mxnest

        if friction:
            self.run_data.geodata.ifriction = 2
        else:
            self.run_data.geodata.ifriction = 0
        
        self.hurricane_data.hurricane_velocity = (velocity * np.cos(angle),velocity * np.sin(angle))
        self.hurricane_data.R_eye_init = eye
        
        self.ml_data.rho = [1025.0,1025.0/rho_ratio]
        self.ml_data.eta = eta

        self.bathy_type = bathy_type

        self.prefix = "ml_angle%s_v%s_m%s_r%s_dof%s" % (
            int(angle * 180.0 / np.pi),int(velocity),mxnest,
            int(rho_ratio*100),int(mx*my))
        
    def write_data_objects(self):
        super(TwoLayerBaseTest,self).write_data_objects()
        
        import topo_data
        topo_data.write_topo_file('topo.data',bathy_type=self.bathy_type,
                                        plot=False,force=False)
        
tests = []
factor = 3

# Two Layer Tests
tests.append(TwoLayerBaseTest(velocity=5.0, angle= 0.00*np.pi, eye=(0.0,0.0), mx=70*factor, my=60*factor))
tests.append(TwoLayerBaseTest(velocity=5.0, angle= 0.25*np.pi, eye=(200e3,-100e3), mx=70*factor, my=60*factor))
tests.append(TwoLayerBaseTest(velocity=5.0, angle=-0.25*np.pi, eye=(200e3,100e3), mx=70*factor, my=60*factor))
tests.append(TwoLayerBaseTest(velocity=5.0, angle= 0.50*np.pi, eye=(425e3,-100e3), mx=70*factor, my=60*factor))
tests.append(TwoLayerBaseTest(velocity=5.0, angle=-0.50*np.pi, eye=(425e3,100e3), mx=70*factor, my=60*factor))

# Speed Tests
for speed in [5.0,10.0,15.0,20.0,30.0]:
    tests.append(TwoLayerBaseTest(velocity=speed, mx=70*factor, my=60*factor))

# R test
# for r in [0.95,0.96,0.97,0.98,0.99]:
#     tests.append(TwoLayerBaseTest(rho_ratio=r, mx=70*factor, my=60*factor))

# Two layer version of shallow sea test
# tests.append(TwoLayerBaseTest(velocity=5.0,angle=0.0,eye=(0.0,0.0),mxnest=1,
#                     mx=70*factor,my=60*factor,bathy_type='flat1000',
#                     friction=False,eta=[0.0,-100.0]))

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'all':
            tests_to_be_run = tests
        else:
            tests_to_be_run = []
            for test in sys.argv[1:]:
                tests_to_be_run.append(tests[int(test)])
            
        test_runs.run_tests(tests_to_be_run,parallel=True)

    else:
        test_runs.print_tests(tests)