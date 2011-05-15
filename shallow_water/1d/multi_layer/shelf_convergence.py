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

import numpy as np
from scipy.linalg import norm
import matplotlib

# matplotlib.rcParams['figure.figsize'] = [7.2,4.0]

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
    
# Parameters
name = 'shelf'
resolutions = [100,200,400,800,1200,1600,2000,4000]
evalue_methods = [1,2,3,4]
index_frame = 150
plot_styles = ['kx','r+','b.']
plot_labels = ['Static Linearized','Dynamic Linearized','Velocity Expansion','LAPACK']
# plot_titles = ["Top Layer Depths","Top Layer Velocities","Bottom Layer Depths","Bottom Layer Velocities"]
base_path = os.path.join(os.environ['DATA_PATH'],"multi_layer_1d",name)

# Error matrix
error = np.zeros((len(resolutions),len(evalue_methods),4))

# Get the LAPACK solutions
solutions_exact = []
run_datas = []
problem_datas = []
for mx in resolutions:
    prefix = "ml_1d_e%s_m%s_%s" % (4,mx,name)
    out_dir = ''.join((prefix,"_output"))
    solutions_exact.append(Solution(index_frame,path=os.path.join(base_path,out_dir)))
    run_datas.append(Data(os.path.join(base_path,out_dir,'claw.data')))
    problem_datas.append(Data(os.path.join(base_path,out_dir,'problem.data')))


# Get the errors
for (m,method) in enumerate(evalue_methods):
    solutions = []
    run_datas = []
    problem_datas = []
    for (i,mx) in enumerate(resolutions):
        # Construct path
        prefix = "ml_1d_e%s_m%s_%s" % (method,mx,name)
        out_dir = ''.join((prefix,"_output"))
    
        sol Solution(index_frame,path=os.path.join(base_path,out_dir))
        run_data = Data(os.path.join(base_path,out_dir,'claw.data'))
        problem_data = Data(os.path.join(base_path,out_dir,'problem.data'))
    
        x,h,u = extract_data(problem_data,sol)
        x,h_exact,u_exact = extract_data(problem_datas[i]],solutions_exact[i])
        
        error[i,m,0] = norm(h[:,0]-h_exact[:,0])
        error[i,m,1] = norm(u[:,0]-u_exact[:,0])
        error[i,m,2] = norm(h[:,1]-h_exact[:,1])
        error[i,m,3] = norm(u[:,1]-u_exact[:,1])

# Top Depth Error Plot
fig = plt.figure()
ax = fig.add_subplot(111)
for (m,method) in enumerate(evalue_methods):
    ax.plot(resolutions,error[:,m,0],plot_style[m])
    
plt.show()

