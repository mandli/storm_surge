#!/usr/bin/env python
# encoding: utf-8
r"""
Contains all test suites for the 1d multi-layer code

Can run a particular test by giving the number of the test at the comand line

:Authors:
    Kyle Mandli (2011-05-11) Initial version
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

tests = []

# Base Plane Wave Test
class PlaneWaveBaseTest(test_runs.TestML2D):
    
    def __init__(self,angle):
        super(PlaneWaveBaseTest,self).__init__()

        self.name = "plane_wave"
        
        self.run_data.clawdata.mx = 100
        self.run_data.clawdata.my = 100
        
        self.ml_data.eigen_method = 2
        self.ml_data.init_type = 2
        self.ml_data.angle = angle

        # Convert angle to degrees for the label
        self.prefix = "ml_2d_angle%s" % int(angle * 180.0 / np.pi)
               
# Angle test suite
for angle in [0.0,np.pi/4.0,np.pi/2.0,-np.pi/4.0]:
    tests.append(PlaneWaveBaseTest(angle))

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'all':
            tests_to_be_run = tests
        else:
            tests_to_be_run = []
            for test in sys.argv[1:]:
                tests_to_be_run.append(tests[int(test)])
            
        test_runs.run_tests(tests,parallel=True)

    else:
        test_runs.print_tests(tests)