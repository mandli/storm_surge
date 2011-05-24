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

class BubbleBaseTest(test_runs.TestML2D):
    
    def __init__(self,eigen_method = 1):
        super(BubbleBaseTest,self).__init__()
        
        self.name = "bubble"
        
        self.run_data.clawdata.mx = 100
        self.run_data.clawdata.my = 100
        
        self.ml_data.eigen_method = eigen_method

# Eigen method tests               
for method in [1,2,3,4]:
    tests.append(BubbleBaseTest(method))

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