#!/usr/bin/env python

import subprocess
import sys

import setrun

# Setup data
run_data = setrun.setrun()
ml_data = setrun.set_multilayer_data()

# Default parameters
run_type = 2
ml_data.inundation_method = 2
ml_data.eigen_method = 2

if len(sys.argv) < 2:
    run_type = 1
    ml_data.inundation_method = 1
elif len(sys.argv) == 2:
    run_type = int(sys.argv[1])
elif len(sys.argv) == 3:
    run_type = int(sys.argv[1])
    ml_data.inundation_method = int(sys.argv[2])
elif len(sys.argv) == 4:
    run_type = int(sys.argv[1])
    ml_data.inundation_method = int(sys.argv[2])
    ml_data.eigen_method = int(sys.argv[3])
    

if run_type == 1:
    run_data.clawdata.outstyle = 3
    run_data.clawdata.iout = [1,100]
    run_data.src_split = 1
    
    ml_data.init_type = 5
    ml_data.init_location = 0.5
    ml_data.eta_1 = 0.0
    ml_data.eta_2 = -0.6
    ml_data.bathy_type = 1
    ml_data.bathy_left = -1.0
    ml_data.bathy_right = -1.0
    ml_data.epsilon = 0.00

    ml_data.wind_type = 0
    
    setplot = 'setplot'

elif run_type == 2:
    run_data.clawdata.mx = 5000
    # run_data.clawdata.outstyle = 3
    # run_data.clawdata.iout = [1,100]
        
    run_data.clawdata.nout = 300
    run_data.clawdata.outstyle = 1
    run_data.clawdata.tfinal = 7200.0
    
    run_data.clawdata.xlower = -400e3
    run_data.clawdata.mthbc_xupper = 3
    
    run_data.src_split = 1
    
    ml_data.rho_air = 1.0
    ml_data.rho_1 = 1025.0
    ml_data.rho_2 = 1028.0
    ml_data.init_type = 4
    ml_data.init_location = 300e3
    ml_data.eta_1 = 0.0
    ml_data.eta_2 = -300.0
    ml_data.epsilon = 0.4
    # ml_data.bathy_location = -30e3
    # ml_data.bathy_left = -4000.0
    # ml_data.bathy_right = -300.0
    ml_data.bathy_type = 2
    ml_data.x0 = -130e3
    ml_data.basin_depth = -4000.0
    ml_data.x1 = -30e3
    ml_data.shelf_depth = -200.0
    ml_data.wind_type = 0
    
    setplot = 'setplot_shelf'
elif run_type == 3:
    
    run_data.clawdata.outstyle = 3
    run_data.clawdata.iout = [1,100]
    run_data.src_split = 1
    
    ml_data.init_type = 6
    ml_data.init_location = 0.5
    ml_data.eta_1 = 0.0
    ml_data.eta_2 = -0.9
    ml_data.bathy_type = 1
    ml_data.bathy_left = -1.0
    ml_data.bathy_right = -1.0
    ml_data.epsilon = 0.25

    ml_data.wind_type = 0
    
    setplot = 'setplot'
    
else:
    print "Invalid run_type =",run_type
    sys.exit(1)

run_data.write()
ml_data.write()

# Run simulation
RUNCLAW_CMD = "python $CLAW/python/pyclaw/runclaw.py"
PLOTCLAW_CMD = "python $CLAW/python/pyclaw/plotters/plotclaw.py"

# output_path = "./%s_%s_%s_output" % (run_type,ml_data.inundation_method,ml_data.eigen_method)
# plots_path = "./%s_%s_%s_plots" % (run_type,ml_data.inundation_method,ml_data.eigen_method)
output_path = "./_output"
plots_path = "./_plots"
data_path = "./"
run_cmd = "%s xclaw %s T F %s" % (RUNCLAW_CMD,output_path,data_path)
plot_cmd = "%s %s %s %s" % (PLOTCLAW_CMD,output_path,plots_path,setplot)
# cmd = ";".join((run_cmd,plot_cmd))
cmd = run_cmd
print cmd
subprocess.Popen(cmd,shell=True).wait()