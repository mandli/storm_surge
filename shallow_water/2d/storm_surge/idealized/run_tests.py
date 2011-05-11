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
    base_path = os.path.join(os.environ['DATA_PATH'],"multi_layer_2d")
else:
    base_path = os.getcwd()
base_path = os.path.expanduser(base_path)
parallel = True
if os.environ.has_key('OMP_NUM_THREADS'):
    max_processes = int(os.environ['OMP_NUM_THREADS'])
else:
    max_processes = 4
process_queue = []
runclaw_cmd = "python $CLAW/python/pyclaw/runclaw.py"
plotclaw_cmd = "python $CLAW/python/pyclaw/plotters/plotclaw.py"
         
test_suites = [{'name':'mx100_3_idealized','setplot':"setplot",
                'run_data':{'mx':100,'my':100},
                'multilayer_data':{'eigen_method':4,'wave_family':3},
                'hurricane_data':{}
               },
               {'name':'mx100_4_idealized','setplot':"setplot",
                'run_data':{'mx':100,'my':100},
                'multilayer_data':{'eigen_method':4,'wave_family':4},
                'hurricane_data':{}
               }]

def run_tests(tests):
    
    for (i,test) in enumerate(tests):
        # Base parameters
        rundata = setrun.setrun()
        h_data = setrun.set_hurricane_data()
        ml_data = setrun.set_multilayer_data()
        
        # Set rundata parameters
        for (key,value) in test['run_data'].iteritems():
            setattr(rundata.clawdata,key,value)
        
        # Set hurricane data
        for (key,value) in test['hurricane_data'].iteritems():
            setattr(h_data,key,value)
        
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
        print "Number of processes currently:",len(process_queue)
        if parallel and len(process_queue) < max_processes - 1:
            process_queue.append(subprocess.Popen(cmd,shell=True,
                stdout=log_file,stderr=log_file))
        else:
            subprocess.Popen(cmd,shell=True,stdout=log_file,
                stderr=log_file).wait()
                
        # Remove any proccess that have completed
        for process in process_queue:
            if process.poll() == 0:
                print "Removing process..."
                process_queue.remove(process)
        print "Number of processes currently:",len(process_queue)
                

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