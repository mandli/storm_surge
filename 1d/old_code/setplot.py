
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

import os

import numpy as np
import re

import matplotlib
import matplotlib.pyplot as mpl

from pyclaw.plotters import geoplot, colormaps

matplotlib.rcParams['figure.figsize'] = [16.0,6.0]

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """

    def hurricane_afterframe(current_data):
        # Draw line for eye of hurricane
        pass

    def bathy(current_data):
        return np.loadtxt(os.path.join(plotdata.outdir,'fort.aux'),converters={0:(lambda x:float(re.compile("[Dd]").sub("e",x)))})
    
    def eta_1(current_data):
        h_1 = current_data.q[:,0]
        return h_1 + bathy(current_data)
        
    def eta_2(current_data):
        h_2 = current_data.q[:,2]
        return h_2 + eta_1(current_data)
        
    def u_1(current_data):
        h_1 = current_data.q[:,0]
        return (np.abs(h_1) > 1e-3) * current_data.q[:,1]/h_1[:]
        
    def u_2(current_data):
        h_2 = current_data.q[:,2]
        return (np.abs(h_2) > 1e-3) * current_data.q[:,3]/h_2[:]

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # ========================================================================
    #  Fill plot
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='full',figno=0)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(1,2,1)'
    plotaxes.title = 'Multilayer Surfaces'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
     
    # Bottom layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta_1
    plotitem.plot_var2 = bathy
    plotitem.color = 'b'
    plotitem.show = True
    
    # Top Layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta_2
    plotitem.plot_var2 = eta_1
    plotitem.color = (0.2,0.8,1.0)
    plotitem.show = True
    
    # Plot bathy
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = bathy
    plotitem.color = 'k'
    # plotitem.plotstyle = '-'
    plotitem.show = True
    
    # Plot line in between layers
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta_1
    plotitem.color = 'k'
    # plotitem.plotstyle = 'o'
    plotitem.show = True
    
    # # Plot line on top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta_2
    plotitem.color = 'k'
    # plotitem.plotstyle = 'x'
    plotitem.show = True
    
    # Layer Velocities
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.title = "Layer Velocities"
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    
    # Bottom layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = u_1
    plotitem.color = 'b'
    plotitem.show = True

    # Top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.color = (0.2,0.8,1.0)
    plotitem.plot_var = u_2
    plotitem.show = True
    
    # # Wind plot
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.axescmd = 'subplot(1,2,2)'
    # plotaxes.title = "Wind Velocity"
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = [-5.0,5.0]
    # 
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = 4
    # plotitem.color = 'r'
    # plotitem.show = True
    

    # ========================================================================
    #  Plot Layer Velocities
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='velocities',figno=1)
    plotfigure.show = True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(1,2,1)'
    plotaxes.title = "Layer Velocities"
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    
    # Bottom layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = u_1
    plotitem.color = 'b'
    plotitem.show = True

    # Top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.color = (0.2,0.8,1.0)
    plotitem.plot_var = u_2
    plotitem.show = True
    
    # Wind plot
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.title = "Wind Velocity"
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 4
    plotitem.color = 'r'
    plotitem.show = True
    
    # ========================================================================
    #  Plot Wind Velocity
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='wind',figno=2)
    plotfigure.show = True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Wind Velocity"
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-5.0,5.0]
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 4
    plotitem.color = 'r'
    plotitem.show = True
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
