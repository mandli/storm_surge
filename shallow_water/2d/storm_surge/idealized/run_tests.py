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

class IdealizedBaseTest(test_runs.TestML2D):
    
    def __init__(self,wave_family,eigen_method):
        super(IdealizedBaseTest,self).__init__()
        
        self.name = "idealized%s" % wave_family
        
        self.run_data.clawdata.mx = 100
        self.run_data.clawdata.my = 100
        self.run_data.clawdata.xlower = 0.0
        self.run_data.clawdata.xupper = 1.0
        self.run_data.clawdata.outstyle = 1
        self.run_data.clawdata.nout = 50
        self.run_data.clawdata.tfinal = 0.5
        
        self.ml_data.eigen_method = eigen_method
        self.ml_data.eta = [0.0,-0.6]
        self.ml_data.rho = [0.95,1.0]
        self.ml_data.wave_family = wave_family
        if 2 <= wave_family <= 3:
            self.ml_data.epsilon = 0.1
        else:
            self.ml_data.epsilon = 0.04

# Method test
for family in [3,4]:
    for method in [1,2,3,4]:
        tests.append(IdealizedBaseTest(family,method))

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