#!/usr/bin/env python
# encoding: utf-8
r"""
Compare 2d idealized runs with 1d

:Authors:
    Kyle Mandli (2011-05-15) Initial version
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
import matplotlib.pyplot as plt

from pyclaw.solution import Solution
from pyclaw.data import Data

def extract_fields_2d(problem_data,sol,slice_index):
    x = sol.grids[0].dimensions[0]
    y = sol.grids[0].dimensions[1]
    h = np.empty((x.n,2))
    u = np.zeros((x.n,2))
    
    h[:,0] = sol.q[:,slice_index,0] / problem_data.rho[0]
    h[:,1] = sol.q[:,slice_index,3] / problem_data.rho[1]
    index = np.nonzero(h[:,0] > 1e-3)
    u[index,0] = sol.q[index,slice_index,1] / sol.q[index,slice_index,0]
    index = np.nonzero(h[:,1] > 1e-3)
    u[index,1] = sol.q[index,slice_index,4] / sol.q[index,slice_index,3]
    
    return x,h,u
    
def extract_fields_1d(problem_data,sol):
    x = sol.grids[0].dimensions[0]
    h = np.empty((x.n,2))
    u = np.zeros((x.n,2))
    
    h[:,0] = sol.q[:,0] / problem_data.rho_1
    h[:,1] = sol.q[:,2] / problem_data.rho_2
    index = np.nonzero(h[:,0] > 1e-3)
    u[index,0] = sol.q[index,1] / sol.q[index,0]
    index = np.nonzero(h[:,1] > 1e-3)
    u[index,1] = sol.q[index,3] / sol.q[index,2]
    
    return x,h,u


# Parameters
plot_styles = ['bx','b+','b.','r.']
plot_labels = ['Static Linearized','Dynamic Linearized','Velocity Expansion','LAPACK']
plot_titles = ["Top Layer Depths","Top Layer Velocities","Bottom Layer Depths","Bottom Layer Velocities"]

frame_index = 25
slice_index = 50
loc = [0,0]

# 3rd wave family comparison
fig = plt.figure()
for field in [0,1,2,3]:
    ax = fig.add_subplot(2,2,field+1)
    
    base_path = os.path.join(os.environ['DATA_PATH'])
    prefix = "ml_%sd_e%s_m%s_%s" % (1,4,500,"idealized_3")
    out_dir = ''.join((prefix,"_output"))
    path_1d = os.path.join(base_path,'multi_layer_1d','idealized_3',out_dir)
        
    sol_1d = Solution(frame_index,path=path_1d)
    prob_data_1d = Data(os.path.join(path_1d,'problem.data'))
    x1,h1,u1 = extract_fields_1d(prob_data_1d,sol_1d)
        
    if field == 0:
        ax.plot(x1.center,h1[:,0],'k-',label='1d Reference')
    elif field == 1:
        ax.plot(x1.center,u1[:,0],'k-',label='1d Reference')
    elif field == 2:
        ax.plot(x1.center,h1[:,1],'k-',label='1d Reference')
    elif field == 3:
        ax.plot(x1.center,u1[:,1],'k-',label='1d Reference')
    
    for method in [1,2,4]:
        prefix = "ml_%sd_e%s_m%s_%s" % (2,method,100,"idealized_redux_3")
        out_dir = ''.join((prefix,"_output"))
        path_2d = os.path.join(base_path,'multi_layer_2d','idealized_redux_3',out_dir)
        
        sol_2d = Solution(frame_index,path=path_2d)
        prob_data_2d = Data(os.path.join(path_2d,'multilayer.data'))
        x2,h2,u2 = extract_fields_2d(prob_data_2d,sol_2d,slice_index)
        
        if field == 0:
            ax.plot(x2.center,h2[:,0],plot_styles[method-1],label=plot_labels[method-1])
        elif field == 1:
            ax.plot(x2.center,u2[:,0],plot_styles[method-1],label=plot_labels[method-1])
        elif field == 2:
            ax.plot(x2.center,h2[:,1],plot_styles[method-1],label=plot_labels[method-1])
        elif field == 3:
            ax.plot(x2.center,u2[:,1],plot_styles[method-1],label=plot_labels[method-1])
    
    ax.set_title(plot_titles[field])
    if field == 0:
        ax.legend(loc=loc[0])

fig.savefig(os.path.join(base_path,'multi_layer_2d','idealized_redux_3','3_comparison.pdf'))

# 4th family comparison
fig = plt.figure()
for field in [0,1,2,3]:
    ax = fig.add_subplot(2,2,field+1)
    
    base_path = os.path.join(os.environ['DATA_PATH'])
    prefix = "ml_%sd_e%s_m%s_%s" % (1,4,500,"idealized_4")
    out_dir = ''.join((prefix,"_output"))
    path_1d = os.path.join(base_path,'multi_layer_1d','idealized_4',out_dir)
        
    sol_1d = Solution(frame_index,path=path_1d)
    prob_data_1d = Data(os.path.join(path_1d,'problem.data'))
    x1,h1,u1 = extract_fields_1d(prob_data_1d,sol_1d)
        
    if field == 0:
        ax.plot(x1.center,h1[:,0],'k-',label='1d Reference')
    elif field == 1:
        ax.plot(x1.center,u1[:,0],'k-',label='1d Reference')
    elif field == 2:
        ax.plot(x1.center,h1[:,1],'k-',label='1d Reference')
    elif field == 3:
        ax.plot(x1.center,u1[:,1],'k-',label='1d Reference')
    
    for method in [1,2,4]:
        prefix = "ml_%sd_e%s_m%s_%s" % (2,method,100,"idealized_redux_4")
        out_dir = ''.join((prefix,"_output"))
        path_2d = os.path.join(base_path,'multi_layer_2d','idealized_redux_4',out_dir)
        
        sol_2d = Solution(frame_index,path=path_2d)
        prob_data_2d = Data(os.path.join(path_2d,'multilayer.data'))
        x2,h2,u2 = extract_fields_2d(prob_data_2d,sol_2d,slice_index)
        
        if field == 0:
            ax.plot(x2.center,h2[:,0],plot_styles[method-1],label=plot_labels[method-1])
        elif field == 1:
            ax.plot(x2.center,u2[:,0],plot_styles[method-1],label=plot_labels[method-1])
        elif field == 2:
            ax.plot(x2.center,h2[:,1],plot_styles[method-1],label=plot_labels[method-1])
        elif field == 3:
            ax.plot(x2.center,u2[:,1],plot_styles[method-1],label=plot_labels[method-1])
    
    ax.set_title(plot_titles[field])

    if field == 0:
        ax.legend(loc=loc[1])

fig.savefig(os.path.join(base_path,'multi_layer_2d','idealized_redux_4','4_comparison.pdf'))