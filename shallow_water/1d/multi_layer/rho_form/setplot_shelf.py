
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

g = 9.81

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """

    def add_dashes(current_data):
        from pylab import ylim,plot
        plot([-30000,-30000], [-1,1],'k--')
        
    def kinetic_energy(current_data):
        q = current_data.q
        h = eta_1(current_data) - bathy(current_data)
        u = 0.5 * (u_1(current_data) + u_2(current_data))
        return 0.5*h*u**2
        
    def potential_energy(current_data):
        return 0.5*eta_1(current_data)**2

    def total_energy(current_data):
        return kinetic_energy(current_data) + potential_energy(current_data)

    def km_labels(current_data):
        r"""Flips xaxis and labels with km"""
        mpl.xlabel('km')
        locs,labels = mpl.xticks()
        labels = np.flipud(locs)/1.e3
        mpl.xticks(locs,labels)
        jump_afteraxes(current_data)

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

    def print_energy(current_data):
        PE = np.sum(potential_energy(current_data))
        KE = np.sum(kinetic_energy(current_data))
        total = PE + KE
        print 'PE = %g, KE = %g, total = %23.16e' % (PE,KE,total)
        
    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.afterframe = print_energy
    
    # ========================================================================
    #  Figure for eta
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='eta', figno=0)
    plotfigure.kwargs = {'figsize':(8,3)}
    plotfigure.show = True

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
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.afteraxes = add_dashes
    #plotaxes.xlimits = [-150.e3, 50e3]
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Energy'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = total_energy
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = [0,1,2]          # list of frames to print
    plotdata.print_fignos = [0,1]            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
