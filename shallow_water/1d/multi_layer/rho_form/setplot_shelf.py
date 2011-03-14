
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

# matplotlib.rcParams['figure.figsize'] = [6.0,10.0]

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """

    def km_labels(current_data):
        r"""Flips xaxis and labels with km"""
        mpl.xlabel('km')
        locs,labels = mpl.xticks()
        labels = np.flipud(locs)/1.e3
        mpl.xticks(locs,labels)
        jump_afteraxes(current_data)

    def jump_afteraxes(current_data):
        # Plot position of jump on plot
        mpl.hold(True)
        mpl.plot([370e3,370e3],[-6,6],'k--')
        mpl.plot([0.0,400e3],[0.0,0.0],'k--')
        mpl.hold(False)

    def bathy(current_data):
        out_dir = current_data.plotdata.outdir
        return np.loadtxt(os.path.join(out_dir,'fort.aux'),converters={0:(lambda x:float(re.compile("[Dd]").sub("e",x)))})
    
    def eta_1(current_data):
        r"""Top surface"""
        h_1 = current_data.q[:,0]
        return h_1 + eta_2(current_data)
        
    def eta_2(current_data):
        r"""Middle surface"""
        h_2 = current_data.q[:,2]
        return h_2 + bathy(current_data)
        
    def u_1(current_data):
        h_1 = current_data.q[:,0]
        index = np.nonzero(h_1 > 1e-3)
        u_1 = np.zeros(h_1.shape)
        u_1[index] = current_data.q[index,1]/h_1[index]
        return u_1
        
    def u_2(current_data):
        h_2 = current_data.q[:,2]
        index = np.nonzero(h_2 > 1e-3)
        u_2 = np.zeros(h_2.shape)
        u_2[index] = current_data.q[index,3] / h_2[index]
        return u_2

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    # Settings
    xlimits = [0.0,400e3]
    ylimits_depth = [-4000.0,100.0]
    xlimits_zoomed = xlimits
    ylimits_surface_zoomed = [-0.5,0.5]
    ylimits_internal_zoomed = [-306,-294]
    ylimits_velocities = [-1.0,1.0]
    
    # ========================================================================
    #  Fill plot
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='full',figno=0)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.title = 'Multilayer Surfaces'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits_depth
    plotaxes.afteraxes = km_labels
     
    # Top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta_1
    plotitem.plot_var2 = eta_2
    plotitem.color = (0.2,0.8,1.0)
    plotitem.show = True
    
    # Bottom Layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta_2
    plotitem.plot_var2 = bathy
    plotitem.color = 'b'
    plotitem.show = True
    
    # Plot bathy
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = bathy
    plotitem.color = 'k'
    # plotitem.plotstyle = '-'
    plotitem.show = True
    
    # Plot line in between layers
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta_2
    plotitem.color = 'k'
    # plotitem.plotstyle = 'o'
    plotitem.show = True
    
    # Plot line on top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta_1
    plotitem.color = 'k'
    # plotitem.plotstyle = 'x'
    plotitem.show = True
    
    # Layer Velocities
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.title = "Layer Velocities"
    plotaxes.xlimits = xlimits
    plotaxes.afteraxes = km_labels
    
    # Bottom layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = u_2
    plotitem.color = 'b'
    plotitem.show = True

    # Top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.color = (0.2,0.8,1.0)
    plotitem.plot_var = u_1
    plotitem.show = True
    
    # ========================================================================
    #  Fill plot zoom
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='full_zoom',figno=1)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,1)'
    plotaxes.title = 'Top Surface'
    plotaxes.xlimits = xlimits_zoomed
    plotaxes.ylimits = ylimits_surface_zoomed
    plotaxes.afteraxes = km_labels
     
    # Top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta_1
    plotitem.plot_var2 = eta_2
    plotitem.color = (0.2,0.8,1.0)
    plotitem.show = True
    
    # Bottom Layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta_2
    plotitem.plot_var2 = bathy
    plotitem.color = 'b'
    plotitem.show = True
    
    # Plot bathy
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = bathy
    plotitem.color = 'k'
    # plotitem.plotstyle = '-'
    plotitem.show = True
    
    # Plot line in between layers
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta_2
    plotitem.color = 'k'
    plotitem.plotstyle = '+'
    plotitem.show = True
    
    # Plot line on top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta_1
    plotitem.color = 'k'
    plotitem.plotstyle = 'x'
    plotitem.show = True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,2)'
    plotaxes.title = 'Internal Surface'
    plotaxes.xlimits = xlimits_zoomed
    plotaxes.ylimits = ylimits_internal_zoomed
    plotaxes.afteraxes = km_labels
     
    # Top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta_1
    plotitem.plot_var2 = eta_2
    plotitem.color = (0.2,0.8,1.0)
    plotitem.show = True
    
    # Bottom Layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_fill_between')
    plotitem.plot_var = eta_2
    plotitem.plot_var2 = bathy
    plotitem.color = 'b'
    plotitem.show = True
    
    # Plot bathy
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = bathy
    plotitem.color = 'k'
    # plotitem.plotstyle = '-'
    plotitem.show = True
    
    # Plot line in between layers
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta_2
    plotitem.color = 'k'
    plotitem.plotstyle = '+'
    plotitem.show = True
    
    # Plot line on top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = eta_1
    plotitem.color = 'k'
    plotitem.plotstyle = 'x'
    plotitem.show = True
    
    # Layer Velocities
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(3,1,3)'
    plotaxes.title = "Layer Velocities"
    plotaxes.xlimits = xlimits_zoomed
    # plotaxes.ylimits = ylimits_velocities
    plotaxes.afteraxes = jump_afteraxes
    
    # Bottom layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = u_2
    plotitem.color = 'b'
    plotitem.plotstyle = '+-'
    plotitem.show = True

    # Top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.color = (0.2,0.8,1.0)
    plotitem.plot_var = u_1
    plotitem.plotstyle = 'x-'
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
    #  h-values
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='depths',figno=2)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.title = 'Depths'
    plotaxes.xlimits = xlimits
    plotaxes.afteraxes = jump_afteraxes
     
    # Top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = (0.2,0.8,1.0)
    plotitem.show = True
    
    # Bottom layer
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 2
    plotitem.color = 'b'
    plotitem.plotstyle = '-'
    plotitem.show = True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.title = 'Depths Zoomed'
    plotaxes.afteraxes = jump_afteraxes
    # plotaxes.xlimits = [0.0,1.0]
    # plotaxes.ylimits = [-1.0,0.5]
    plotaxes.xlimits = xlimits_zoomed
    # plotaxes.xlimits = [0.0,2000.0]
    # plotaxes.ylimits = [-2000.0,100.0]
     
    # Top layer
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'x'
    plotitem.color = (0.2,0.8,1.0)
    plotitem.show = True
    
    # Bottom layer
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 2
    plotitem.color = 'b'
    plotitem.plotstyle = '+'
    plotitem.show = True
    
    # ========================================================================
    #  Plot Wind Velocity
    # ========================================================================
    # plotfigure = plotdata.new_plotfigure(name='wind',figno=2)
    # plotfigure.show = True
    # 
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.title = "Wind Velocity"
    # plotaxes.xlimits = 'auto'
    # plotaxes.ylimits = [-5.0,5.0]
    # 
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = 4
    # plotitem.color = 'r'
    # plotitem.show = True
    
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

    
