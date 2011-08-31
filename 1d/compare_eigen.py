#!/usr/bin/env python
# encoding: utf-8
r"""
Compare the last results assuming they are different eigenvalue methods

:Authors:
    Kyle Mandli (2011-05-14) Initial version
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
import subprocess

import numpy as np
from scipy.linalg import norm
import matplotlib

factor = 0.8
matplotlib.rcParams['figure.figsize'] = [7.2,4.0]

import matplotlib.pyplot as plt

from pyclaw.data import Data
from pyclaw.solution import Solution

def extract_data(problem_data,sol):
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
    
def reformat_title(title):
    result = title.split()[0].lower()

    for word in title.split()[1:]:
        result = "_".join((result,word.lower()))
        
    return result


# Parameters
name = 'shelf'
mx = 2000
index_frame = 150
plot_styles = ['bx','b+','b.','k-']
plot_labels = ['Static Linearized','Dynamic Linearized','Velocity Expansion','LAPACK']
plot_titles = ["Top Layer Depths","Top Layer Velocities","Bottom Layer Depths","Bottom Layer Velocities"]
base_path = os.path.join(os.environ['DATA_PATH'],"multi_layer_1d",name)
solutions = []
run_datas = []
problem_datas = []

# Load data
for method in [1,2,3,4]:
    # Construct path
    prefix = "ml_1d_e%s_m%s_%s" % (method,mx,name)
    out_dir = ''.join((prefix,"_output"))
    
    solutions.append(Solution(index_frame,path=os.path.join(base_path,out_dir)))
    run_datas.append(Data(os.path.join(base_path,out_dir,'claw.data')))
    problem_datas.append(Data(os.path.join(base_path,out_dir,'problem.data')))


# Setup exact solution for error estimates
x_exact,h_exact,u_exact = extract_data(problem_datas[2],solutions[2])

# Plots
figs = []
axes = []
for i in [0,1,2,3]:
    figs.append(plt.figure(i))
    axes.append(figs[i].add_subplot(111))

table_string = ""
for (i,sol) in enumerate(solutions):
    x,h,u = extract_data(problem_datas[i],sol)
    
    print "Method: ",plot_labels[i]
    print "Top layer - depth difference       =",norm(h[:,0] - h_exact[:,0])
    print "Top layer - velocity difference    =",norm(u[:,0] - u_exact[:,0])
    print "Bottom layer - depth difference    =",norm(h[:,1] - h_exact[:,1])
    print "Bottom layer - velocity difference =",norm(u[:,1] - u_exact[:,1])
    error_string = "%s & %s & %s & %s & %s" % (i+1,norm(h[:,0] - h_exact[:,0]),
                                                   norm(u[:,0] - u_exact[:,0]),
                                                   norm(h[:,1] - h_exact[:,1]),
                                                   norm(u[:,1] - u_exact[:,1]))
    table_string = "\\\\ \n".join((table_string,error_string))
    
    # Top layer - depth
    axes[0].plot(x.center,h[:,0],plot_styles[i],label=plot_labels[i])
    # Top layer - velocity
    axes[1].plot(x.center,u[:,0],plot_styles[i],label=plot_labels[i])
    # Bottom layer - depth
    axes[2].plot(x.center,h[:,1],plot_styles[i],label=plot_labels[i])
    # Bottom layer - velocity
    axes[3].plot(x.center,u[:,1],plot_styles[i],label=plot_labels[i])
    
print table_string
    
# ===============   =============
# Location String   Location Code
# ===============   =============
# 'best'            0
# 'upper right'     1
# 'upper left'      2
# 'lower left'      3
# 'lower right'     4
# 'right'           5
# 'center left'     6
# 'center right'    7
# 'lower center'    8
# 'upper center'    9
# 'center'          10
# ===============   =============

axes[0].legend(loc=3)

for i in [0,1,2,3]:
    axes[i].set_title(plot_titles[i])
    figs[i].savefig(os.path.join(base_path,"%s.pdf" % reformat_title(plot_titles[i])))

# Also plot the contours
if name == "shelf" and True:
    mx = 4000
    for method in [1,2,3,4]:
        prefix = "ml_1d_e%s_m%s_%s" % (method,mx,name)
        out_dir = os.path.join(base_path,''.join((prefix,"_output")))
        cmd = "python shelf_contour.py %s" % out_dir
        subprocess.Popen(cmd,shell=True).wait()
