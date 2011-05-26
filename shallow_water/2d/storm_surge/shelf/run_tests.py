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

import sys
import os

import numpy as np

import test_runs
         
tests = []

class ShelfBaseTest(test_runs.TestML2D):
    
    def __init__(self,eigen_method):
        super(ShelfBaseTest,self).__init__()
        
        self.name = "shelf"
                
        self.run_data.clawdata.xlower = -400000.0
        self.run_data.clawdata.xupper = 0.0
        self.run_data.clawdata.ylower = -300e3
        self.run_data.clawdata.yupper = 300e3
        self.run_data.clawdata.mx = 500
        self.run_data.clawdata.my = 60
        
        self.run_data.clawdata.nout = 300
        self.run_data.clawdata.tfinal = 7200.0
        
        self.run_data.clawdata.mthbc_xupper = 3
        
        self.run_data.clawdata.mxnest = -1
        self.run_data.geodata.topofiles = []
        
        self.ml_data.eigen_method = eigen_method
        self.ml_data.rho = [1025.0,1028.0]
        
        self.ml_data.init_type = 5
        self.ml_data.init_location = [300e3,50e3]
        self.ml_data.wave_family = wave_family
        self.ml_data.epsilon = 0.4
        self.ml_data.sigma = 0.02
        
        self.ml_data.bathy_type = 1
        self.ml_data.bathy_location = -30e3
        self.ml_data.bathy_left = -4000
        self.ml_data.bathy_right = -200
        
        self.prefix = "ml_e%s" % eigen_method

# Eigenmethod tests
for method in [1,2,3,4]:
    tests.append(ShelfBaseTest(method))


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