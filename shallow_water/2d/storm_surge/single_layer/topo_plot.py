#!/usr/bin/env python
# encoding: utf-8
r"""
Plot profile and overview plot of topography

:Authors:
    Kyle Mandli (2010-07-09) Initial version
"""

import os

import numpy as np
import matplotlib.pyplot as plt

import pyclaw.util as util
import pyclaw.plotters.geoplot as geoplot
import pyclaw.plotters.gaugetools as gt
import pyclaw.plotters.colormaps as colormaps

# Generate bathy
N = 100
x = np.linspace(-200e3,500e3,N)
y = np.linspace(-300e3,300e3,N)
[X,Y] = np.meshgrid(x,y)
shallow_shelf = [(477e3,-100),(500e3,80.0)]
topo = util.create_topo_func(shallow_shelf)
bath = np.empty(X.shape)
for i in xrange(len(x)):
    for j in xrange(len(y)):
        bath[i,j] = topo(x[i],y[j])

# Profile plot
topo_file = '/Users/mandli/Documents/research/Papers/awr10/figures/ss_topo.png'
if not os.path.exists(topo_file):
    plt.figure(1)
    index = np.trunc(N/2.0)
    plt.hold(True)
    plt.plot(x,bath[:,index],'k',linewidth=2)
    plt.fill_between(x,bath[:,index],x*0.0,where=(bath[:,index] < 0.0),color=[0,0.5,1])
    plt.axis([-200e3,500e3,-120,100])
    plt.title("Bathymetry Cross-Section")
    plt.xlabel('km')
    plt.ylabel('m')
    locs,labels = plt.xticks()
    labels = locs/1e3
    plt.xticks(locs,labels)
    plt.hold(False)
    plt.savefig(topo_file)
    plt.show()

# Domain plot
domain_file = '/Users/mandli/Documents/research/Papers/awr10/figures/ss_domain.png'
# domain_file = './ss_domain.png'
if not os.path.exists(domain_file):
    plt.figure(2)
    plt.hold(True)
    
    # Water and land
    water = np.ma.masked_where(bath>0.0,np.zeros(bath.shape))
    plt.pcolor(X,Y,bath.T,cmap=geoplot.land_colors)
    # light blue = [0,0.5,1]
    plt.pcolor(X,Y,water.T,cmap=colormaps.make_colormap({1.0:'w',0.0:'w'}))
    
    # Hurricane Track
    plt.plot([-200e3,500e3],[0.0,0.0],'r-',lw=2)
    for n in [0,5]:
        plt.plot(n*3600*5,0.0,'rD',markersize=10)
        plt.text(n*3600*5-4e3,15e3,"%s hr" % n)
    for n in [-6,13,20]:
        plt.plot(n*3600*5,0.0,'rD',markersize=10)
        plt.text(n*3600*5-10e3,15e3,"%s hr" % n)
    # plt.annotate('Position of Storm (h)',(0.0,-5e3),(-100e3,-100e3),arrowprops={'width':1,'color':'k','shrink':0.10})

    # Gauges
    setgauges = gt.read_setgauges('./_output')
    offset = 30e3
    n = 9
    xn = setgauges.x[n]
    yn = setgauges.y[n]
    plt.plot([xn],[yn],'ko')
    plt.text(xn-offset,yn-8e3,str(1))
    n = 11
    xn = setgauges.x[n]
    yn = setgauges.y[n]
    plt.plot([xn],[yn],'ko')
    plt.text(xn-offset,yn-8e3,str(2))
    # n = 18
    # xn = setgauges.x[n]
    # yn = setgauges.y[n]
    # plt.plot([xn],[yn],'ko')
    # plt.text(xn-offset,yn-8e3,str(3))
    
    # General settings
    plt.axis([-200e3,500e3,-300e3,300e3])
    plt.title("Shallow Sea Domain")
    plt.xlabel('km')
    plt.ylabel('km')
    locs,labels = plt.xticks()
    labels = locs/1e3
    plt.xticks(locs,labels)
    locs,labels = plt.yticks()
    labels = locs/1e3
    plt.yticks(locs,labels)
    plt.hold(False)
    plt.savefig(domain_file)
    plt.show()