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
import copy
import glob

import numpy as np

from pyclaw.runclaw import runclaw
from pyclaw.plotters.plotclaw import plotclaw

import setrun

# Parameters
if os.environ.has_key('DATA_PATH'):
    base_path = os.path.join(os.environ['DATA_PATH'],"multi_layer_2d")
else:
    base_path = os.getcwd()
base_path = os.path.expanduser(base_path)
parallel = True
poll_interval = 15.0
if os.environ.has_key('OMP_NUM_THREADS'):
    max_processes = int(os.environ['OMP_NUM_THREADS'])
else:
    max_processes = 4
process_queue = []
runclaw_cmd = "python $CLAW/python/pyclaw/runclaw.py"
plotclaw_cmd = "python $CLAW/python/pyclaw/plotters/plotclaw.py"
         
test_suites = []

base_shelf_test = {'name':'shelf_redux','setplot':'setplot',
                   'run_data':{'mx':2000,'my':100,'nout':300,'outstyle':1,
                      'tfinal':7200.0,'xlower':-400000.0,'mthbc_xupper':3},
                   'multilayer_data':{'rho_air':1.0,'rho_1':1025.0,'rho_2':1028.0,
                      'eigen_method':1,'init_type':5,'init_location':300e3,
                      'eta_2':-300,'epsilon':0.4,'bathy_location':-30e3,
                      'bathy_left':-4000,'bathy_right':-200,'wind_type':0}
                  }

# Eigenmethod tests
for method in [1,2,3,4]:
    test = copy.deepcopy(base_shelf_test)
    test['multilayer_data']['eigen_method'] = method
    test_suites.append(test)

def run_tests(tests):
    
    for (i,test) in enumerate(tests):
        # Create base parameters
        rundata = setrun.setrun()
        ml_data = setrun.set_multilayer_data()
        hurricane_data = setrun.set_hurricane_data()
        
        # Set rundata parameters
        for (key,value) in test['run_data'].iteritems():
            setattr(rundata.clawdata,key,value)
        
        # Set multilayer data
        for (key,value) in test['multilayer_data'].iteritems():
            setattr(ml_data,key,value)
            
        for (key,value) in test['hurricane_data'].iteritems():
            setattr(hurricane_data,key,value)
        
        # Create output paths
        prefix = "ml_%sd_e%s_m%s_%s" % (rundata.clawdata.ndim,
                                        ml_data.eigen_method,
                                        rundata.clawdata.mx,
                                        test['name'])
        data_dirname = ''.join((prefix,'_data'))
        output_dirname = ''.join((prefix,"_output"))
        plots_dirname = ''.join((prefix,"_plots"))
        log_name = ''.join((prefix,"_log"))

        data_path = os.path.join(base_path,test['name'],data_dirname)
        output_path = os.path.join(base_path,test['name'],output_dirname)
        plots_path = os.path.join(base_path,test['name'],plots_dirname)
        log_path = os.path.join(base_path,test['name'],log_name)
        
        # Create test directory if not present
        if not os.path.exists(os.path.join(base_path,test['name'])):
            os.mkdir(os.path.join(base_path,test['name']))
        
        # Clobber old data directory
        if os.path.exists(data_path):
            data_files = glob.glob(os.path.join(data_path,'*.data'))
            for data_file in data_files:
                os.remove(data_file)
        else:
            os.mkdir(data_path)
        
        # Open and start log file
        log_file = open(log_path,'w')
        tm = time.localtime()
        year = str(tm[0]).zfill(4)
        month = str(tm[1]).zfill(2)
        day = str(tm[2]).zfill(2)
        hour = str(tm[3]).zfill(2)
        minute = str(tm[4]).zfill(2)
        second = str(tm[5]).zfill(2)
        date = '%s/%s/%s - %s:%s.%s\n' % (year,month,day,hour,minute,second)
        log_file.write(date)
        
        # Write out data files to output directory
        temp_location = os.getcwd()
        os.chdir(data_path)
        rundata.write()
        ml_data.write()
        hurricane_data.write()
        os.chdir(temp_location)

        # Run the simulation
        run_cmd = "%s xclaw %s T F %s" % (runclaw_cmd,output_path,data_path)
        plot_cmd = "%s %s %s %s" % (plotclaw_cmd,output_path,plots_path,test['setplot'])
        tar_cmd = "tar -cvzf %s.tgz %s" % (plots_path,plots_path)
        cmd = ";".join((run_cmd,plot_cmd))
        # cmd = run_cmd
        print cmd
        if parallel:
            print "Number of processes currently:",len(process_queue)
            while len(process_queue) == max_processes:
                for process in process_queue:
                    if process.poll() == 0:
                        process_queue.remove(process)
                time.sleep(poll_interval)
            process_queue.append(subprocess.Popen(cmd,shell=True,
                stdout=log_file,stderr=log_file))
            
        else:
            # subprocess.Popen(cmd,shell=True).wait()
            subprocess.Popen(cmd,shell=True,stdout=log_file,
                stderr=log_file).wait()
                

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
        if sys.argv[1].lower() == 'all':
            tests = test_suites
        else:
            tests = []
            for test in sys.argv[1:]:
                tests.append(test_suites[int(test)])
        run_tests(tests)
    else:
        print_tests()