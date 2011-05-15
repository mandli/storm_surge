
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
from pyclaw.data import Data

# matplotlib.rcParams['figure.figsize'] = [6.0,10.0]

def setplot(plotdata):
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """
    claw_data = Data(os.path.join(plotdata.outdir,'claw.data'))
    prob_data = Data(os.path.join(plotdata.outdir,'problem.data'))
    g = 9.81
        
    # ========================================================================
    #  Plot variable functions
    def bathy(current_data):
        out_dir = current_data.plotdata.outdir
        return np.loadtxt(os.path.join(out_dir,'fort.aux'),
            converters={0:(lambda x:float(re.compile("[Dd]").sub("e",x)))})
    
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

    def kinetic_energy(current_data):
        q = current_data.q
        h = eta_1(current_data) - bathy(current_data)
        u = 0.5 * (u_1(current_data) + u_2(current_data))
        return 0.5*h*u**2
        
    def potential_energy(current_data):
        return 0.5*eta_1(current_data)**2

    def total_energy(current_data):
        return kinetic_energy(current_data) + potential_energy(current_data)

    def print_energy(current_data):
        PE = np.sum(potential_energy(current_data))
        KE = np.sum(kinetic_energy(current_data))
        total = PE + KE
        print 'PE = %g, KE = %g, total = %23.16e' % (PE,KE,total)
        
            
    # ========================================================================
    #  Labels    
    def add_bathy_dashes(current_data):
        mpl.hold(True)
        mpl.plot([-30e3,-30e3],[-10,10],'k--')
        mpl.hold(False)
        
    def add_horizontal_dashes(current_data):
        mpl.hold(True)
        mpl.plot([-400e3,0.0],[0.0,0.0],'k--')
        mpl.hold(False)

    def km_labels(current_data):
        r"""Flips xaxis and labels with km"""
        mpl.xlabel('km')
        locs,labels = mpl.xticks()
        labels = np.flipud(locs)/1.e3
        mpl.xticks(locs,labels)
        
    def time_labels(current_data):
        r"""Convert time to hours"""
        pass
        
    
    # ========================================================================
    # Limit Settings
    xlimits = [claw_data.xlower,claw_data.xupper]
    ylimits_depth = [-4000.0,100.0]
    xlimits_zoomed = [-30e3-1e3,-30e3+1e3]
    ylimits_surface_zoomed = [prob_data.eta_1 - 0.5,prob_data.eta_1 + 0.5]
    ylimits_internal_zoomed = [prob_data.eta_2 - 2.5,prob_data.eta_2 + 2.5] 
    # ylimits_velocities = [-1.0,1.0]
    ylimits_velocities = [-0.04,0.04]
    ylimits_kappa = [0.0,1.2]
        
    # Create data object
    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.afterframe = print_energy
    
    # ========================================================================
    #  Figure for eta
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='eta', figno=0)
    plotfigure.kwargs = {'figsize':(8,3)}
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([0.1,0.15,0.80,0.75])'
    def fixfig(current_data):
        from pylab import xticks,yticks,xlabel,ylabel,savefig,ylim,title
        t = current_data.t
        add_dashes(current_data)
        xticks([-300000,-200000,-100000, -30000],['300','200','100','30','0'],\
          fontsize=15)
        ylim(-0.4,0.6)
        yticks([-0.4,-0.2,0,0.2,0.4],fontsize=15)
        #xlabel('kilometres offshore', fontsize=15)
        ylabel('Metres', fontsize=15)
        title('Surface at t = %i seconds' % int(t),fontsize=20)
        # savefig('shelf%s.eps' % str(current_data.frameno).zfill(2))

    plotaxes.afteraxes = fixfig
    #plotaxes.xlimits = [-150.e3, 50e3]
    #plotaxes.ylimits = [-0.4, 0.4]
    plotaxes.title = 'Surface'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    def surface(current_data):
        eta = B + current_data.q[:,0]
        return eta
    plotitem.plot_var = eta_1
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':2}
    plotitem.show = True       # show on plot?
    
    # ========================================================================
    #  Figure for energy
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='energy', figno=1)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Energy'
    def energy_axes(current_data):
        add_horizontal_dashes(current_data)
        add_bathy_dashes(current_data)
        km_labels(current_data)
    plotaxes.afteraxes = energy_axes

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = total_energy
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotitem.show = True 
    
    # ========================================================================
    #  Fill plot zoom
    # ========================================================================
    def fill_items(plotaxes):
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
        plotitem.show = True
            
        # Plot line in between layers
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.plot_var = eta_2
        plotitem.color = 'k'
        plotitem.plotstyle = '-'
        plotitem.show = True
    
        # Plot line on top layer
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.plot_var = eta_1
        plotitem.color = 'k'
        plotitem.plotstyle = '-'
        plotitem.show = True
    
    # Top surface
    plotfigure = plotdata.new_plotfigure(name='full_zoom',figno=100)
    plotfigure.show = False
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Top Surface'
    plotaxes.xlimits = xlimits_zoomed
    plotaxes.ylimits = ylimits_surface_zoomed
    def top_afteraxes(current_data):
        km_labels(current_data)
        add_bathy_dashes(current_data)
    plotaxes.afteraxes = top_afteraxes
     
    plotaxes = fill_items(plotaxes)

    # Internal surface
    plotfigure = plotdata.new_plotfigure(name='internal_zoom',figno=101)
    plotfigure.show = False
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Internal Surface'
    plotaxes.xlimits = xlimits_zoomed
    plotaxes.ylimits = ylimits_internal_zoomed
    plotaxes.afteraxes = km_labels
     
    plotaxes = fill_items(plotaxes)
    
    # ========================================================================
    #  Velocities
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name="Velocities",figno=200)
    plotfigure.show = False
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Layer Velocities"
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits_velocities
    def velocity_afteraxes(cd):
        add_bathy_dashes(cd)
        add_horizontal_dashes(cd)
        # km_labels(cd)
        mpl.title("Layer Velocities t = %4.1f s" % cd.t)
        mpl.xticks([-300e3,-200e3,-100e3,-30e3],[300,200,100,30],fontsize=15)
        mpl.xlabel('km')
    plotaxes.afteraxes = velocity_afteraxes
    
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
    
    # ========================================================================
    #  Velocities with Kappa
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='vel_kappa',figno=14)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize':(7,6)}
    
    def twin_axes(cd):
        fig = mpl.gcf()
        fig.clf()
        
        x = cd.grid.dimensions[0].center
        
        # Draw velocity and kappa plot
        ax1 = fig.add_subplot(111)     # the velocity scale
        ax2 = ax1.twinx()              # the kappa scale
        
        # Bottom layer velocity
        bottom_layer = ax1.plot(x,u_2(cd),'k-',label="Bottom Layer Velocity")
        # Top Layer velocity
        top_layer = ax1.plot(x,u_1(cd),'b-',label="Top Layer velocity")#,color=(0.2,0.8,1.0))
        
        # Kappa
        kappa_line = ax2.plot(x,cd.q[:,5],color='r',label="Kappa")
        ax2.plot(x,np.ones(x.shape),'r--')

        ax1.set_xlabel('km')
        mpl.xticks([-300e3,-200e3,-100e3,-30e3],[300,200,100,30],fontsize=15)

        ax1.plot([prob_data.bathy_location,prob_data.bathy_location],ylimits_velocities,'k--')
        ax1.legend((bottom_layer,top_layer,kappa_line),('Bottom Layer','Top Layer',"Kappa"),loc=3)
        ax1.set_title("Layer Velocities and Kappa t = %4.1f s" % cd.t)
        ax1.set_ylabel('Velocities (m/s)')
        ax2.set_ylabel('Kappa (1/Ri)')
        ax1.set_xlim(xlimits)
        ax1.set_ylim(ylimits_velocities)
        ax2.set_ylim(ylimits_kappa)
        
        # mpl.subplots_adjust(hspace=0.1)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.afteraxes = twin_axes
    
    # ========================================================================
    #  Combined Top and Internal Surface
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='combined_surface',figno=13)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize':(6,6)}
    
    # Top surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.title = 'Surfaces'
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits_surface_zoomed
    def top_afteraxes(cd):
        mpl.xlabel('')
        locs,labels = mpl.xticks()
        # labels = np.flipud(locs)/1.e3
        labels = ['' for i in xrange(len(locs))]
        mpl.xticks(locs,labels)
        add_bathy_dashes(cd)
        mpl.ylabel('m')
        mpl.title("Surfaces t = %4.1f s" % cd.t)
    plotaxes.afteraxes = top_afteraxes
    plotaxes = fill_items(plotaxes)
    
    # Internal surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.title = ''
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits_internal_zoomed
    def internal_surf_afteraxes(cd):
        km_labels(cd)
        mpl.title('')
        mpl.ylabel('m')
        mpl.subplots_adjust(hspace=0.05)
        mpl.xticks([-300e3,-200e3,-100e3,-30e3],[300,200,100,30],fontsize=15)
        mpl.xlabel('km')
    plotaxes.afteraxes = internal_surf_afteraxes
    plotaxes = fill_items(plotaxes)
    
    
    
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

    
