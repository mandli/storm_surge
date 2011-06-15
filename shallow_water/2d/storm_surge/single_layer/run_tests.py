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

class SingleLayerBaseTest(test_runs.TestML2D):
    
    def __init__(self,velocity=5.0,angle=0.0,eye=(0.0,0.0),mxnest=5,
                    mx=70,my=60,bathy_type='simple_shelf',friction=True):
        super(SingleLayerBaseTest,self).__init__()
        
        self.type = "storm_surge"
        self.name = "single_layer"
        self.setplot = "setplot"
        
        self.run_data.clawdata.mx = mx
        self.run_data.clawdata.my = my
        self.run_data.clawdata.mxnest = -mxnest
        self.hurricane_data.hurricane_velocity = (velocity * np.cos(angle),velocity * np.sin(angle))
        self.hurricane_data.R_eye_init = eye

        if friction:
            self.run_data.geodata.ifriction = 2
        else:
            self.run_data.geodata.ifriction = 0

        self.bathy_type = bathy_type

        self.prefix = "sl_angle%s_m%s_v%s_dof%s" % (int(angle * 180.0 / np.pi),mxnest,int(velocity),int(mx*my))
        

    
    def write_data_objects(self):
        super(SingleLayerBaseTest,self).write_data_objects()
        
        import topo_data
        topo_data.write_topo_file('topo.data',bathy_type=self.bathy_type,
                                        plot=False,force=False)

tests = []

# Angle Tests
# tests.append(SingleLayerBaseTest(5.0,0.0,(0.0,0.0),False))
# tests.append(SingleLayerBaseTest(5.0,0.25*np.pi,(200e3,-100e3),False))
# tests.append(SingleLayerBaseTest(5.0,-0.25*np.pi,(200e3,100e3),False))
# tests.append(SingleLayerBaseTest(5.0,-0.50*np.pi,(425e3,100e3),False))
# tests.append(SingleLayerBaseTest(5.0, 0.50*np.pi,(425e3,-100e3),False))

# Comparison test case
levels = 3
tests.append(SingleLayerBaseTest(velocity=5.0,angle=0.0,eye=(0.0,0.0),mxnest=1,
                    mx=70*levels,my=60*levels,bathy_type='shallow_shelf',
                    friction=False))

# Speed Tests
for speed in [5.0,10.0,15.0,20.0,30.0]:
    tests.append(SingleLayerBaseTest(speed,0.0,(0.0,0.0),True))

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