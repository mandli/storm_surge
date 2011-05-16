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

matplotlib.rcParams['figure.figsize'] = [6.0,4.0]

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
resolutions = [100,200,400,800,1200,1600,2000,3000,4000,5000]
evalue_methods = [1,2,3]
index_frame = 67
legend_location = [0,0,0,0]
plot_style = ['ko-','rs-','bd-']
plot_label = ['Static Linearized','Dynamic Linearized','Velocity Expansion','LAPACK']
plot_title = ["Top Layer Depths","Top Layer Velocities","Bottom Layer Depths","Bottom Layer Velocities"]
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
    for (i,mx) in enumerate(resolutions):
        # Construct path
        prefix = "ml_1d_e%s_m%s_%s" % (method,mx,name)
        out_dir = ''.join((prefix,"_output"))
    
        sol = Solution(index_frame,path=os.path.join(base_path,out_dir))
        run_data = Data(os.path.join(base_path,out_dir,'claw.data'))
        problem_data = Data(os.path.join(base_path,out_dir,'problem.data'))
    
        x,h,u = extract_data(problem_data,sol)
        x,h_exact,u_exact = extract_data(problem_datas[i],solutions_exact[i])
        
        error[i,m,0] = norm(h[:,0]-h_exact[:,0])
        error[i,m,1] = norm(u[:,0]-u_exact[:,0])
        error[i,m,2] = norm(h[:,1]-h_exact[:,1])
        error[i,m,3] = norm(u[:,1]-u_exact[:,1])

# Top Depth Error Plot
figures = []
for field in [0,1,2,3]:
    figures.append(plt.figure())
    ax = figures[field].add_subplot(1,1,1)
    for (m,method) in enumerate(evalue_methods):
        ax.plot(resolutions,error[:,m,field],plot_style[m],label=plot_label[m])
    
    ax.set_title(plot_title[field])
    ax.legend(loc=legend_location[field])
    ax.set_xlabel('N')
    ax.set_ylabel("Error")
    locs,labels = plt.xticks()
    labels = [100,800,1600,4000]
    locs = labels
    plt.xticks(locs,labels)
    
    fname = "%s_convergence.pdf" % reformat_title(plot_title[field])
    figures[field].savefig(os.path.join(base_path,fname))

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

# Construct tables
for (m,method) in enumerate(evalue_methods):
    table_section = ""
    for (i,mx) in enumerate(resolutions):
        if i == 0:
            table_row = "\\\\ \n\multirow{8}{*}{%s} & %s" % (m+1,mx)
        else:
            table_row = "\\\\ \n & %s" % mx
        for field in [0,1,2,3]:
            table_row = " & ".join((table_row,str(error[i,m,field])))
        table_section = "".join((table_section,table_row))
    print "Method: %s" % plot_label[m]
    print table_section
