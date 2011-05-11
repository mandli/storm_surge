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
import time

import numpy as np

from pyclaw.runclaw import runclaw
from pyclaw.plotters.plotclaw import plotclaw

import setrun

# Parameters
if os.environ.has_key('DATA_PATH'):
    base_path = os.path.join(os.environ['DATA_PATH'],"multi_layer_1d")
else:
    base_path = os.getcwd()
base_path = os.path.expanduser(base_path)
parallel = False
process_queue = []
runclaw_cmd = "python $CLAW/python/pyclaw/runclaw.py"
plotclaw_cmd = "python $CLAW/python/pyclaw/plotters/plotclaw.py"
         
test_suites = [{'name':'rp_waves_3','setplot':'setplot',
                'run_data':{'xlower':0.0,'xupper':1.0,'mx':500,'outstyle':1,
                    'tfinal':1.0},
                'multilayer_data':{'eigen_method':1,'init_type':1,
                    "init_location":0.45,"wave_family":3,'eta_2':-0.6,
                    'bathy_left':-1.0,'bathy_right':-0.2,'wind_type':0}
                },
                {'name':'rp_waves_4','setplot':'setplot',
                'run_data':{'xlower':0.0,'xupper':1.0,'mx':500,'outstyle':1,
                    'tfinal':0.1},
                'multilayer_data':{'eigen_method':1,'init_type':1,
                    "init_location":0.45,"wave_family":3,'eta_2':-0.6,
                    'bathy_left':-1.0,'bathy_right':-0.2,'wind_type':0}
                },
                {'name':'oscillatory_wind','setplot':'setplot_oscillatory',
                 'run_data':{'mx':100,'outstyle':1,'nout':160,'tfinal':10.0,
                    'mthbc_xlower':3,'mthbc_xupper':3},
                 'multilayer_data':{'init_type':0,'eta_1':0.0,'eta_2':-0.25,
                    'wind_type':3,'A':5.0,"rho_air":1.15,"rho_1":1025,
                    "rho_2":1045,"N":2.0,'omega':2.0,"t_length":10.0,
                    'bathy_left':-1.0,'bathy_right':-1.0}
                },
                {'name':'shelf','setplot':'setplot_shelf',
                 'run_data':{'mx':2000,'nout':300,'outstyle':1,'tfinal':7200.0,
                    'xlower':-400000.0,'mthbc_xupper':3},
                 'multilayer_data':{'rho_air':1.0,'rho_1':1025.0,'rho_2':1028.0,
                    'eigen_method':4,'init_type':4,'init_location':300e3,
                    'eta_2':-300,'epsilon':0.4,'bathy_location':-30e3,
                    'bathy_left':-4000,'bathy_right':-200,'wind_type':0}
                }]

def run_tests(tests):
    
    for (i,test) in enumerate(tests):
        # Base parameters
        rundata = setrun.setrun()
        ml_data = setrun.set_multilayer(rundata)
        
        # Set rundata parameters
        for (key,value) in test['run_data'].iteritems():
            setattr(rundata.clawdata,key,value)
        
        # Set multilayer data
        for (key,value) in test['multilayer_data'].iteritems():
            setattr(ml_data,key,value)
        
        # Write out data files
        rundata.write()
        ml_data.write()
        
        # Create output directory
        prefix = "ml_1d_e%s_%s" % (ml_data.eigen_method,test['name'])
        tm = time.localtime()
        year = str(tm[0]).zfill(4)
        month = str(tm[1]).zfill(2)
        day = str(tm[2]).zfill(2)
        hour = str(tm[3]).zfill(2)
        minute = str(tm[4]).zfill(2)
        second = str(tm[5]).zfill(2)
        date = ''
        #date = '_%s%s%s-%s%s%s' % (year,month,day,hour,minute,second)
        output_dirname = ''.join((prefix,date,"_output"))
        plots_dirname = ''.join((prefix,date,"_plots"))
        log_name = ''.join((prefix,date,"_log"))

        output_path = os.path.join(base_path,output_dirname)
        plots_path = os.path.join(base_path,plots_dirname)
        log_path = os.path.join(base_path,log_name)

        log_file = open(log_path,'w')

        # Run the simulation
        run_cmd = "%s xclaw %s" % (runclaw_cmd,output_path)
        plot_cmd = "%s %s %s %s" % (plotclaw_cmd,output_path,plots_path,test['setplot'])
        tar_cmd = "tar -cvzf %s.tgz %s" % (plots_path,plots_path)
        cmd = ";".join((run_cmd,plot_cmd))
        print cmd
        if parallel:
            process_queue.append(subprocess.Popen(cmd,shell=True,
                stdout=log_file,stderr=log_file))
        else:
            subprocess.Popen(cmd,shell=True).wait()

def print_tests():
    for (i,test) in enumerate(test_suites):
        print "Test %s: %s" % (i,test['name'])
        print "  Setplot: %s" % test['setplot']
        print "  Run Data:"
        print "    %s" % test['run_data']
        print "  Multilayer Data:"
        print "    %s" % test['multilayer_data']

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1].lowercase() == 'all':
            tests = test_suites
        else:
            tests = []
            for test in sys.argv[1:]:
                tests.append(test_suites[int(test)])
        run_tests(tests)
    else:
        print_tests()