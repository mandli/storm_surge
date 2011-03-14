
import os
import re
import numpy as np
import pylab

from pyclaw.plotters.data import ClawPlotData

# Plot settings
pd = ClawPlotData()
pd.outdir = "_output"

# Load the bathymetry
b =  np.loadtxt(os.path.join(pd.outdir,'fort.aux'),
                converters={0:(lambda x:float(re.compile("[Dd]").sub("e",x)))})

def read_data():
    num_frames = 30
    mx = 2000
    eta = np.ndarray((mx,num_frames,2))
    t = np.ndarray((num_frames))
    
    for frameno in xrange(num_frames):
        frame = pd.getframe(frameno)
        q = frame.grids[0].q
        t[frameno] = frame.t / 3600.0
        eta[:,frameno,1] = q[:,2] + b
        eta[:,frameno,0] = q[:,0] + eta[:,frameno,1]
    
    x = frame.p_center[0]
    X,T = np.meshgrid(x,t)
    return X,T,eta


def contour_plot(X,T,eta):
    pylab.figure(1,figsize=[7,8])
    pylab.clf()
    pylab.axes([.1,.1,.6,.8])
    clines = np.linspace(.025,.4,15)
    pylab.contour(X,T,eta[:,:,0].T,clines,colors='r')
    # pylab.contour(X,T,eta,-clines,colors='b',linestyles='solid')
    pylab.contour(X,T,eta[:,:,0].T,-clines,colors='b')
    pylab.plot([-30e3,-30e3],[0,2],'k--')
    pylab.xticks([-300e3,-200e3,-100e3,-30e3],[300,200,100,30],fontsize=15)
    pylab.yticks(fontsize=15)
    pylab.xlabel("Kilometers offshore",fontsize=20)
    pylab.ylabel("Hours",fontsize=20)
    pylab.title("Contours of surface",fontsize=20)
    add_timeslices()
    pylab.savefig("shelf_top.png")
    
    # pylab.figure(2,figsize=[7,8])
    # pylab.clf()
    # pylab.axes([.1,.1,.6,.8])
    # clines = np.linspace(.025,.4,15)
    # pylab.contour(X,T,eta[:,:,1].T,clines,colors='r')
    # #pylab.contour(X,T,eta,-clines,colors='b',linestyles='solid')
    # pylab.contour(X,T,eta[:,:,1].T,-clines,colors='b')
    # pylab.plot([-30e3,-30e3],[0,2],'k--')
    # pylab.xticks([-300e3,-200e3,-100e3,-30e3],[300,200,100,30],fontsize=15)
    # pylab.yticks(fontsize=15)
    # pylab.xlabel("Kilometres offshore",fontsize=20)
    # pylab.ylabel("Hours",fontsize=20)
    # pylab.title("Contours of surface",fontsize=20)
    # add_timeslices()
    # pylab.savefig("shelf_internal.png")


def add_timeslices():
    times = pylab.array([0,200,400,600,1000,1400,2000,2800,3400,4800]) 
    for t in times:
        thours = t / 3600.
        pylab.plot([-400e3,0.],[thours,thours],'k')
        pylab.annotate('%s seconds' % t,[0.,thours],[30e3,thours],\
              arrowprops={'width':1,'color':'k','frac':0.2,'shrink':0.1})


def xt_and_frameplots(X,T,eta):
    """
    Unfinished attempt to plot frames with xt plot
    """
    contour_plot(X,T,eta)
    framenos = [0,1,2,3,5,7,10,14,27,24]
    plotdata = ClawPlotData()
    plotdata = setplot(plotdata)
    for frameno in framenos:
        frametools.plotframe(frameno,plotdata)

def mesh_plot(X,T,eta):
    from enthought.mayavi import mlab
    X = 2000*X/np.abs(x.min())
    eta = 1000*eta
    mlab.figure(1)
    mlab.clf()
    mlab.mesh(X,T,eta)

if __name__=="__main__":
    X,T,eta = read_data()
    contour_plot(X,T,eta)
    # pylab.savefig("shelfxt.eps")
