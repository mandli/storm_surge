#!/usr/bin/env python

import sys
import os
import subprocess
import glob

figures = [0,1]

# Figure parameters
delay = 10
loop = 1

# Expand figure paths
if len(sys.argv) > 1:
    figure_paths = [os.path.abspath(path) for path in sys.argv[1:]]
else:
    figure_paths = os.path.abspath("./")

movie_dir = os.path.expandvars("$DATA_PATH/storm_surge/movies")

for figure_path in figure_paths:
    sim_name = os.path.split(figure_path)[-1][:-6]
    for figure in figures:
        movie_name = "".join((sim_name,"_fig",str(figure),".gif"))
        plot_glob = os.path.join(figure_path,"frame*fig%s.png" % figure)
        cmd="convert -delay %s -loop %s %s %s" % (delay,loop,plot_glob,os.path.join(movie_dir,movie_name))
        print(cmd)
        subprocess.Popen(cmd,shell=True).wait()
