#!/usr/bin/env python

import os

import numpy as np
import matplotlib.pyplot as plt

import compare_gauges

# Base directory
base_path = os.path.expandvars("$DATA_PATH/storm_surge/")
output_path = os.path.join(base_path,"gauge_comparisons")
if not os.path.exists(output_path):
    os.mkdir(output_path)

# Comparisons
comparisons = [["single_layer/sl_angle0_v5_m5_dof4200_output","multi_layer/ml_angle0_v5_m1_r97_dof37800_output"],
               ["single_layer/sl_angle0_v10_m5_dof4200_output","multi_layer/ml_angle0_v10_m1_r97_dof37800_output"],
               ["single_layer/sl_angle0_v15_m5_dof4200_output","multi_layer/ml_angle0_v15_m1_r97_dof37800_output"],
               ["single_layer/sl_angle0_v20_m5_dof4200_output","multi_layer/ml_angle0_v20_m1_r97_dof37800_output"],
               ["single_layer/sl_angle0_v30_m5_dof4200_output","multi_layer/ml_angle0_v30_m1_r97_dof37800_output"],
               ["single_layer/sl_angle-45_v5_m5_dof4200_output","multi_layer/ml_angle-45_v5_m1_r97_dof37800_output"],
               ["single_layer/sl_angle45_v5_m5_dof4200_output","multi_layer/ml_angle45_v5_m1_r97_dof37800_output"],
               ["single_layer/sl_angle-90_v5_m5_dof4200_output","multi_layer/ml_angle-90_v5_m1_r97_dof37800_output"],
               ["single_layer/sl_angle90_v5_m5_dof4200_output","multi_layer/ml_angle90_v5_m1_r97_dof37800_output"]]

comparison_paths = []
for path in comparisons:
    comparison_paths.append([os.path.join(base_path,path[0]),os.path.join(base_path,path[1])])
              
for (i,data_paths) in enumerate(comparison_paths):
    figures = compare_gauges.plot_gauges('all',data_paths,plot_vars=[[3],[6]],
                               titles=["Top Surface","Bottom Surface"],
                               kwargs={'figsize':(8,4)})
    for (num,figure) in figures.iteritems():
        name = os.path.split(data_paths[0])[-1][:-7]
        figure.savefig(os.path.join(output_path,"%s_%s_gauge.pdf" % (name,num)))

figure = compare_gauges.plot_gauge_locations(
            os.path.join(base_path,"single_layer/sl_angle0_v5_m5_dof4200_output"),
            offset=[40e3,10e3],bathy_ref_lines=[350e3,450e3,480e3])
            
ax = figure.axes[0]
ax.set_title("Gauge Locations")
ax.set_xticks([-200e3,-100e3,0,100e3,200e3,300e3,400e3,500e3])
ax.set_yticks([-200e3,-100e3,0,100e3,200e3])
ax.set_xticklabels([-200,-100,0,100,200,300,400,500])
ax.set_yticklabels([-200,-100,0,100,200])
ax.set_xlabel('km')
ax.set_ylabel('km')
plt.grid(True)

figure.savefig(os.path.join(output_path,"gauge_locations.pdf"))
