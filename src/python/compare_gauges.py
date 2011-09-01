#!/usr/bin/env python
r"""Plot the requested gauges verse each other"""

import sys
import os

import numpy as np
import matplotlib.pyplot as plt

import pyclaw.data

class TideGauge(object):
    
    def __init__(self,number=None,times=None,data=None,location=None):
        self.number = number
        self.t = times
        self.q = data
        self.location = location
    
    def plot_gauge(self,ax,meqn=0,style='b-'):
        ax.plot(self.t,self.q[:,meqn],style)
    
    def plot_gauge_location(self,ax,style='ko'):
        ax.plot(self.location[0],self.location[1],style)
        
    def __str__(self):
        output =  "Gauge %s: \n" % self.number
        output += "  Location = (%s,%s)\n" % (self.location[0],self.location[1])
        output += "  t = [%s,%s]\n" % (self.t[0],self.t[-1])
        output += "  q.shape = %s\n" % str(self.q.shape)
        return output
        
    def __repr__(self):
        return str(self)
        
        

def read_setgauges(path):
    gauge_file = open(os.path.join(path,'setgauges.data'),'r')
    # Read through header
    for i in xrange(6):
        gauge_file.readline()
    num_gauges = int(gauge_file.readline().split()[0])
    # Gauge locations
    gauges = np.empty((num_gauges,3))
    for i in xrange(num_gauges):
        data = gauge_file.readline().split()
        gauges[i,:] = data[0:3]
        
    return gauges

def plot_gauges(gauge_numbers,paths,plot_vars=[[0]],titles=["q(0)"],
                tlim=None,kwargs={}):
    
    if isinstance(paths,basestring):
        paths = [paths]
    
    figures = {}
    styles = ['b','r.','kx','g+']
    
    # Read in each file
    for (k,path) in enumerate(paths):
        print "Reading from %s" % path
        # Read setgauge.data file
        gauge_locations = read_setgauges(path)
        
        gauges = []
        if gauge_numbers == 'all':
            gauge_numbers = [int(num) for num in gauge_locations[:,0]]
            
        # Create gauges
        data = np.loadtxt(os.path.join(path,'fort.gauge'))
        for (i,num) in enumerate(gauge_numbers):
            print "Parsing gauge %s" % num
            # Add location and number
            index = gauge_locations[:,0] == num
            location = gauge_locations[index]
            gauges.append(TideGauge(number=int(location[0,0]),location=location[0,1:]))
            
            # Extract gauge data
            index = data[:,0] == num
            gauge_data = data[index,...]
            gauges[i].t = np.array(gauge_data[:,2])
            gauges[i].q = np.array(gauge_data[:,3:])
        
            # Create a figure for this gauge number unless it already exists
            if num not in figures.keys():
                figures[num] = plt.figure(**kwargs)
                for j in xrange(len(plot_vars[k])):
                    figures[num].add_subplot(1,len(plot_vars[k]),j+1)
                    title = "%s - Gauge %s" % (titles[j],num) 
                    figures[num].axes[j].set_title(title)
                    figures[num].axes[j].set_xlim(tlim)
                
            # Plot the gauge on the appropriate figure
            for (j,m) in enumerate(plot_vars[k]):
                gauges[i].plot_gauge(figures[num].axes[j],meqn=m,style=styles[k])
                
    # for (num,figure) in figures.iteritems():
    #     figure.savefig("./%s_%s" % (num,base_output_name))
    
    return figures
            
def plot_gauge_locations(path,file_name='amr2ez.data',offset=[0,0],bathy_ref_lines=[]):
    # Open data file to find domain
    run_data = pyclaw.data.Data(os.path.join(path,file_name))
    
    # Extract domain parameters
    mx = run_data.mx
    my = run_data.my
    x = np.linspace(run_data.xlower,run_data.xupper,mx)
    y = np.linspace(run_data.ylower,run_data.yupper,my)

    # Read in gauge data
    gauges = read_setgauges(path)
    
    # Plot gauges
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in xrange(gauges.shape[0]):
        ax.plot(gauges[i,1],gauges[i,2],'ko',markersize=5)
        ax.text(gauges[i,1]-offset[0],gauges[i,2]-offset[1],str(int(gauges[i,0])).zfill(2),fontsize=15)
        
    # Plot bathy locations
    for bathy_ref in bathy_ref_lines:
        ax.plot(bathy_ref*np.ones(y.shape),y,'k--')
    
    ax.set_xbound((run_data.xlower,run_data.xupper))
    ax.set_ybound((run_data.ylower,run_data.yupper))
    
    return fig

if __name__ == "__main__":
    if len(sys.argv) > 2:
        gauge_numbers = sys.argv[1]
        data_paths = []
        for path in sys.argv[2:]:
            data_paths.append(os.path.abspath(path))
        plot_gauges(gauge_numbers,data_paths)
    else:
        plot_gauges('all',os.path.abspath('./'))
        