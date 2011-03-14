
from pyclaw.plotters.data import ClawPlotData
import numpy as np
import pylab


def read_data():
    pd = ClawPlotData()
    pd.outdir = "_output_xt"
    fname = pd.outdir + '/fort.H'
    B = np.loadtxt(fname)
    print "Loaded B"

    eta = []
    times = []
    for frameno in range(301):
        frame = pd.getframe(frameno)
        q = frame.grids[0].q
        t = frame.t
        eta.append(B+q[:,0])
        times.append(t/3600.)
    
    x = frame.p_center[0]
    x = x
    X,T = np.meshgrid(x,times)
    eta = np.array(eta)
    return X,T,eta

def contour_plot(X,T,eta):
    pylab.figure(1,figsize=[7,8])
    pylab.clf()
    pylab.axes([.1,.1,.6,.8])
    clines = np.linspace(.025,.4,15)
    pylab.contour(X,T,eta,clines,colors='r')
    #pylab.contour(X,T,eta,-clines,colors='b',linestyles='solid')
    pylab.contour(X,T,eta,-clines,colors='b')
    pylab.plot([-30e3,-30e3],[0,2],'k--')
    pylab.xticks([-300e3,-200e3,-100e3,-30e3],[300,200,100,30],fontsize=15)
    pylab.yticks(fontsize=15)
    pylab.xlabel("Kilometres offshore",fontsize=20)
    pylab.ylabel("Hours",fontsize=20)
    pylab.title("Contours of surface",fontsize=20)

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
    add_timeslices()
    pylab.savefig("shelfxt.png")
    pylab.savefig("shelfxt.eps")
