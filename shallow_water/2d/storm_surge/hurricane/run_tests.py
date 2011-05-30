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
    
    def __init__(self,velocity=5.0,angle=0.0,eye=(0.0,0.0),mxnest=5):
        super(TwoLayerBaseTest,self).__init__()
        
        # Create topography
        
        self.type = "storm_surge"
        self.name = "multi_layer"
        self.setplot = "setplot"
        
        self.run_data.clawdata.mxnest = -mxnest
        self.hurricane_velocity = (velocity * np.cos(angle),velocity * np.sin(angle))
        self.R_eye_init = eye
        
        self.prefix = "ml_angle%s_m%s" % (int(angle * 180.0 / np.pi),mxnest)
        
    def write_data_objects(self):
        super(TwoLayerBaseTest,self).write_data_objects()
        
        import topo_data
        topo_data.write_topo_file('topo.data',bathy_type='simple_shelf',
                                        plot=False,force=False)
        
tests = []

# Two Layer Tests
tests.append(TwoLayerBaseTest(5.0,0.0,(0.0,0.0),5))
tests.append(TwoLayerBaseTest(5.0,0.25*np.pi,(200e3,100e3),5))
tests.append(TwoLayerBaseTest(5.0,-0.25*np.pi,(200e3,-100e3),5))
tests.append(TwoLayerBaseTest(5.0,-0.50*np.pi,(425e3,100e3),5))
tests.append(TwoLayerBaseTest(5.0, 0.50*np.pi,(425e3,-100e3),5))
tests.append(TwoLayerBaseTest(0.0,0.0,(0.0,0.0),5))

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