#!/usr/bin/env python
# encoding: utf-8
r"""
Utilities for running clawpack tests

:Authors:
    Kyle Mandli (2011-05-18) Initial version
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
import copy
import time
import glob

import numpy as np

from pyclaw.runclaw import runclaw
from pyclaw.plotters.plotclaw import plotclaw

def print_tests(tests):
    for (i,test) in enumerate(tests):
        print "====== Test #%s ===========================" % i
        print str(test)

def run_tests(tests,data_path=None,max_processes=None,parallel=True,
                verbose=False,terminal_output=False):
    # Parameters
    if data_path is None:
        if os.environ.has_key('DATA_PATH'):
            base_path = os.path.join(os.environ['DATA_PATH'],test[0].name)
        else:
            base_path = os.getcwd()
    base_path = os.path.expanduser(base_path)
    poll_interval = 15.0
    if max_processes is None:
        if os.environ.has_key('OMP_NUM_THREADS'):
            max_processes = int(os.environ['OMP_NUM_THREADS'])
        else:
            max_processes = 4
    
    process_queue = []
    RUNCLAW_CMD = "python $CLAW/python/pyclaw/runclaw.py"
    PLOTCLAW_CMD = "python $CLAW/python/pyclaw/plotters/plotclaw.py"
    
    # Run tests
    for (i,test) in enumerate(tests):
        # Write out data files locally
        test.create_data_objects()
        test.set_data_objects()
    
        # Create output directory
        data_dirname = ''.join((prefix,'_data'))
        output_dirname = ''.join((test.prefix,"_output"))
        plots_dirname = ''.join((test.prefix,"_plots"))
        log_name = ''.join((test.prefix,"_log.txt"))

        output_path = os.path.join(base_path,output_dirname)
        plots_path = os.path.join(base_path,plots_dirname)
        log_path = os.path.join(base_path,log_name)

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
        date = 'Started %s/%s/%s-%s:%s.%s' % (year,month,day,hour,minute,second)
        log_file.write(date)
    
        # Write out data
        temp_path = os.getcwd()
        os.chdir(data_path)
        test.write_data_objects(data_objects)
        os.chdir(temp_path)
        
        # Run the simulation
        run_cmd = "%s xclaw %s T F %s" % (runclaw_cmd,output_path,data_path)
        plot_cmd = "%s %s %s %s" % (plotclaw_cmd,output_path,plots_path,test['setplot'])
        tar_cmd = "tar -cvzf %s.tgz %s" % (plots_path,plots_path)
        cmd = ";".join((run_cmd,plot_cmd))
        if verbose:
            print cmd
        if parallel:
            while len(process_queue) == max_processes:
                if verbose:
                    print "Number of processes currently:",len(process_queue)
                for process in process_queue:
                    if process.poll() == 0:
                        process_queue.remove(process)
                time.sleep(poll_interval)
            process_queue.append(subprocess.Popen(cmd,shell=True,
                    stdout=log_file,stderr=log_file))
            
        else:
            if terminal_output:
                log_file.write("Outputting to terminal...")
                subprocess.Popen(cmd,shell=True).wait()
                log_file.write("Command completed.")
            else:
                subprocess.Popen(cmd,shell=True,stdout=log_file,
                    stderr=log_file).wait()
    
    # Wait to exit while processes are still going
    while len(process_queue) > 0:
        for process in process_queue:
            if process.poll() == 0:
                process_queue.remove(process)
                print "Number of processes currently:",len(process_queue)
        time.sleep(poll_interval)


class Test(object):
    
    def __init__(self):
        # Base test traits
        self.name = "test"
        self.prefix = ""
        self.setplot = "setplot"
        
    def __str__(self):
        output = "Test %s:" % self.name
        output += "\n  Setplot: %s" % self.setplot
        return output
        
    def set_data_objects(self):
        pass
        
    def write_data_objects(self):
        pass

class TestML2D(Test):
    
    def __init__(self):
        super(TestML2D,self).__init__()
        
        import setrun
    
        self.run_data = setrun.setrun()
        self.ml_data = setrun.set_multilayer_data()
        self.hurricane_data = setrun.set_hurricane_data()
        
        self.prefix = "ml_1d_e%s_m%s_%s" % (self.ml_data.eigen_method,
                                            self.run_data.clawdata.mx,
                                            self.name)
                                            
    def __str__(self):
        output = super(TestML2D,self).__str__()
        output += "\nRundata:    \n%s" % self.run_data
        output += "\n--------------------------------------"
        output += "\nMultilayer: \n%s" % self.ml_data
        output += "\n--------------------------------------"
        output += "\nHurricane:  \n%s" % self.hurricane_data
        return output
    
    def set_data_objects(self):
        # Fill in data object differences
        #
        # Reset prefix
        self.prefix = "ml_1d_e%s_m%s_%s" % (self.ml_data.eigen_method,
                                            self.run_data.clawdata.mx,
                                            self.name)
    
    def write_data_objects(self):
        self.run_data.write()
        self.ml_data.write()
        self.hurricane_data.write()



class TestML1D(Test):
    
    def __init__(self):
        super(TestML1D,self).__init__()
        
        import setrun
        
        self.run_data = setrun.setrun()
        self.ml_data = setrun.set_multilayer_data()
        
        self.prefix = "ml_1d_e%s_m%s_%s" % (self.ml_data.eigen_method,
                                            self.run_data.clawdata.mx,
                                            self.name)
                                            
    def __str__(self):
        output = super(TestML1D,self).__str__()
        output += "\n  Rundata:    \n%s" % self.run_data
        output += "\n--------------------------------------"
        output += "\n  Multilayer: \n%s" % self.ml_data
        return output
    
    def set_data_objects(self):
        # Fill in data object differences
        
        # Reset prefix
        self.prefix = "ml_1d_e%s_m%s_%s" % (self.ml_data.eigen_method,
                                            self.run_data.clawdata.mx,
                                            self.name)
    
    def write_data_objects(self):
        self.run_data.write()
        self.ml_data.write()
        
    

        