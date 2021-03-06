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

import subprocess
import sys
import os
import numpy as np

import test_runs

tests = []

class SmoothBaseTest(test_runs.TestML1D):

    def __init__(self,wave_family,mx=500,epsilon=0.1):
        
        super(IdealizedBaseTest,self).__init__()
        
        self.name = "smooth_%s" % wave_family
        
        # Static data
        self.run_data.clawdata.xlower = 0.0
        self.run_data.clawdata.xupper = 1.0
        self.run_data.clawdata.outstyle = 1
        self.run_data.clawdata.nout = 50
        self.run_data.clawdata.tfinal = 0.5
        
        self.ml_data.init_type = 1
        self.ml_data.init_location = 0.45
        self.ml_data.eta = [0.0,-0.6]
        self.ml_data.bathy_type = 2
        self.ml_data.bathy_left = -1.0
        self.ml_data.bathy_right = -1.0
        self.ml_data.wind_type = 0
        self.ml_data.rho = [0.95,1.0]
        
        # Parameters
        self.run_data.clawdata.mx = mx
        self.ml_data.wave_family = wave_family
        self.ml_data.eigen_method = eigen_method
        self.ml_data.epsilon = epsilon
        
        self.prefix = "ml_e%s_m%s" % (eigen_method,mx)

class IdealizedBaseTest(test_runs.TestML1D):
    
    def __init__(self,wave_family,mx=500,eigen_method=2,epsilon=0.1):
        
        super(IdealizedBaseTest,self).__init__()
        
        self.name = "idealized_%s" % wave_family
        
        # Static data
        self.run_data.clawdata.xlower = 0.0
        self.run_data.clawdata.xupper = 1.0
        self.run_data.clawdata.outstyle = 1
        # self.run_data.clawdata.iout = [1,100]
        self.run_data.clawdata.nout = 50
        self.run_data.clawdata.tfinal = 0.5
        
        self.ml_data.init_type = 1
        self.ml_data.init_location = 0.45
        self.ml_data.eta = [0.0,-0.6]
        self.ml_data.bathy_left = -1.0
        self.ml_data.bathy_right = -0.2
        self.ml_data.wind_type = 0
        self.ml_data.rho = [0.95,1.0]
        
        # Parameters
        self.run_data.clawdata.mx = mx
        self.ml_data.wave_family = wave_family
        self.ml_data.eigen_method = eigen_method
        self.ml_data.epsilon = epsilon
        
        self.prefix = "ml_e%s_m%s" % (eigen_method,mx)

        
class OscillatoryWindBaseTest(test_runs.TestML1D):
    
    def __init__(self,mx=100,eigen_method=2):
        
        super(OscillatoryWindBaseTest,self).__init__()
        
        self.name = "oscillatory_wind"
        self.setplot = "setplot_oscillatory"
        
        self.run_data.clawdata.mx = mx
        self.run_data.clawdata.outstyle = 1
        self.run_data.clawdata.nout = 160
        self.run_data.clawdata.tfinal = 10.0
        self.run_data.clawdata.mthbc_xlower = 3
        self.run_data.clawdata.mthbc_xupper = 3
        
        self.ml_data.init_type = 0
        self.ml_data.eta = [0.0,-0.25]
        self.ml_data.wind_type = 3
        self.ml_data.A = 5.0
        self.ml_data.rho_air = 1.15
        self.ml_data.rho = [1025.0,1045.0]
        self.ml_data.N = 2.0
        self.ml_data.omega = 2.0
        self.ml_data.t_length = (self.run_data.clawdata.tfinal 
                                    - self.run_data.clawdata.t0)
        self.ml_data.bathy_left = -1.0
        self.ml_data.bathy_right = -1.0
        self.ml_data.eigen_method = eigen_method
        
        self.prefix = "ml_e%s_m%s" % (eigen_method,mx)
        
        
class ShelfBaseTest(test_runs.TestML1D):
    
    def __init__(self,eigen_method=2,mx=2000):
        
        super(ShelfBaseTest,self).__init__()
        
        # Test parameters
        self.name = "shelf"
        self.setplot = "setplot_shelf"

        # Data parameters
        self.ml_data.eigen_method = eigen_method
        self.run_data.clawdata.mx = mx
        
        self.run_data.clawdata.nout = 300
        self.run_data.clawdata.outstyle = 1
        self.run_data.clawdata.tfinal = 7200.0
        self.run_data.clawdata.xlower = -400e3
        self.run_data.clawdata.mthbc_xupper = 3
        
        
        self.ml_data.inundation_method = 0
        self.ml_data.rho_air = 1.0
        self.ml_data.rho = [1025.0,1028.0]
        self.ml_data.init_type = 4
        self.ml_data.init_location = 300e3
        self.ml_data.eta = [0.0,-300.0]
        self.ml_data.epsilon = 0.4
        self.ml_data.bathy_type = 1
        self.ml_data.bathy_location = -30e3
        self.ml_data.bathy_left = -4000.0
        self.ml_data.bathy_right = -100.0
        self.ml_data.wind_type = 0
        
        self.prefix = "ml_e%s_m%s" % (eigen_method,mx)
        
class RealShelfBaseTest(ShelfBaseTest):
    
    def __init__(self,eigen_method=2,mx=2000):
        
        super(RealShelfBaseTest,self).__init__(eigen_method,mx)
    
        self.name = "real_shelf"
        self.setplot = "setplot_shelf"
        
        self.ml_data.inundation_method = 2
        self.ml_data.bathy_type = 2
        
        self.prefix = "ml_e%s_m%s" % (eigen_method,mx)

class DryStateBaseTest(test_runs.TestML1D):

    def __init__(self,mx=500,eigen_method=2,entropy_fix=True):

        super(DryStateBaseTest,self).__init__()

        self.name = "dry_state"
        self.setplot = "setplot_drystate"

        self.ml_data.eigen_method = eigen_method
        self.run_data.clawdata.mx = mx

        self.run_data.clawdata.iout = [1,100]
        self.run_data.clawdata.outstyle = 3
        self.run_data.clawdata.xlower = 0.0
        self.run_data.clawdata.xupper = 1.0

        self.ml_data.inundation_method = 2
        self.ml_data.manning = 0.0
        self.ml_data.rho = [0.95,1.0]
        self.ml_data.eta_left = [0.0,-0.5]
        self.ml_data.eta_right = [0.0,-1.0]
        self.ml_data.init_location = 0.5
        
        self.ml_data.bathy_type = 1
        self.ml_data.bathy_left = -1.0
        self.ml_data.bathy_right = -1.0
        self.ml_data.bathy_location = 0.5
        self.ml_data.wind_type = 0

        self.ml_data.entropy_fix = entropy_fix

        self.prefix = "ml_e%s_m%s_fix" % (eigen_method,mx)
        if entropy_fix:
            self.prefix = "".join((self.prefix,"T"))
        else:
            self.prefix = "".join((self.prefix,"F"))
            
class RarefactionBaseTest(test_runs.TestML1D):
    
    def __init__(self,mx=500,eigen_method=2,entropy_fix=True):
        
        super(RarefactionBaseTest,self).__init__()
        
        self.name = "all_rare"
        self.setplot = "setplot_drystate"
        
        self.ml_data.eigen_method = eigen_method
        self.run_data.clawdata.mx = mx

        self.run_data.clawdata.iout = [1,100]
        self.run_data.clawdata.outstyle = 3
        self.run_data.clawdata.xlower = 0.0
        self.run_data.clawdata.xupper = 1.0

        self.ml_data.inundation_method = 2
        self.ml_data.manning = 0.0
        self.ml_data.rho = [0.95,1.0]
        self.ml_data.init_type = 0
        self.ml_data.init_location = 0.5
        self.ml_data.eta = [0.0,-0.5]
        self.ml_data.u_left = [0.0,-0.5]
        self.ml_data.u_right = [0.0,0.5]
        self.ml_data.bathy_type = 1
        self.ml_data.bathy_left = -1.0
        self.ml_data.bathy_right = -1.0
        self.ml_data.bathy_location = 0.5
        self.ml_data.wind_type = 0
        
        self.ml_data.entropy_fix = entropy_fix
        
        self.prefix = "ml_e%s_m%s_fix" % (eigen_method,mx)
        if entropy_fix:
            self.prefix = "".join((self.prefix,"T"))
        else:
            self.prefix = "".join((self.prefix,"F"))
    


# Idealized 3 eigen_method test
for method in [1,2,3,4]:
    tests.append(IdealizedBaseTest(3,epsilon=0.1,eigen_method=method))
     
# Idealized 4 eigen_method test
for method in [1,2,3,4]:
    tests.append(IdealizedBaseTest(4,epsilon=0.04,eigen_method=method))

# Eigen method tests for oscillatory wind
for method in [1,2,3,4]:
    tests.append(OscillatoryWindBaseTest(eigen_method=method))
    
# Convergence test for shelf
# for method in [1,2,3,4]:
    # for mx in [100,200,400,800,1200,1600,2000,3000,4000,5000]:

# This shelf jump case does not fully work with inundation solver for all
# eigenspace methods
mx = 2000
for method in [1,2,3,4]:
    tests.append(ShelfBaseTest(mx=mx,eigen_method=method))
    tests.append(RealShelfBaseTest(mx=mx,eigen_method=method))
    
# Entropy correction tests
tests.append(DryStateBaseTest(eigen_method=2,mx=500,entropy_fix=False))
for method in [1,2,3,4]:
    tests.append(DryStateBaseTest(eigen_method=method,mx=500,entropy_fix=True))

for method in [2,4]:
    tests.append(RarefactionBaseTest(eigen_method=method,entropy_fix=False))
    tests.append(RarefactionBaseTest(eigen_method=method,entropy_fix=True))



if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'all':
            tests_to_be_run = tests
        else:
            tests_to_be_run = []
            for test in sys.argv[1:]:
                tests_to_be_run.append(tests[int(test)])
            
        test_runs.run_tests(tests_to_be_run,parallel=True,plot=True)

    else:
        test_runs.print_tests(tests)