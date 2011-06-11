
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

def meters_to_km(x,y):
    return x/1e3,y/1e3
    
#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    import os

    import numpy as np
    import matplotlib.pyplot as plt

    from pyclaw.plotters import colormaps, geoplot
    from pyclaw.data import Data

    amrdata = Data(os.path.join(plotdata.outdir,'amr2ez.data'))
    hurricane_data = Data(os.path.join(plotdata.outdir,'hurricane.data'))
    multilayer_data = Data(os.path.join(plotdata.outdir,'multilayer.data'))
    
    if multilayer_data.bathy_type == 1:
        ref_lines = [multilayer_data.bathy_location]
    elif multilayer_data.bathy_type == 2:
        ref_lines = [multilayer_data.x0,multilayer_data.x1,multilayer_data.x2]
    else:
        ref_lines = []
    
    plotdata.clearfigures()
    plotdata.clear_frames = False
    plotdata.clear_figs = True
    
    plotdata.save_frames = False
    
    # ========================================================================
    #  Generic helper functions
    # ========================================================================        
    def bathy_ref_lines(current_data):
        plt.hold(True)
        y = [amrdata.ylower,amrdata.yupper]
        for ref_line in ref_lines:
            plt.plot([ref_line,ref_line],y,'y--')
        plt.hold(False)
        
    def day_figure_title(current_data):
        t = current_data.t
        title = current_data.plotaxes.title
        plt.title('%s at time t = %s days' % (title,str(t/(3600.0*24.0))))
        
    def m_to_km_labels(current_data=None):
        plt.xlabel('km')
        plt.ylabel('km')
        locs,labels = plt.xticks()
        labels = locs/1.e3
        plt.xticks(locs,labels)
        locs,labels = plt.yticks()
        labels = locs/1.e3
        plt.yticks(locs,labels)
    
    def pcolor_afteraxes(current_data):
        day_figure_title(current_data)
        m_to_km_labels()
        bathy_ref_lines(current_data)
        # gauge_locations(current_data)
        
    def contour_afteraxes(current_data):
        day_figure_title(current_data)
        m_to_km_labels()
        bathy_ref_lines(current_data)

    # ========================================================================
    # Gauge functions
    # ========================================================================
    def gauge_locations(current_data,gaugenos='all'):
        from pyclaw.plotters import gaugetools
        plt.hold(True)
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=gaugenos, format_string='kx', add_labels=True)
        plt.hold(False)

    def gaugetopo(current_data):
        q = current_data.q
        h = q[:,0]
        eta = q[:,3]
        topo = eta - h
        return topo
        
    def gauge_afteraxes(current_data):
        # Change time to hours
        plt.xlabel('t (days)')
        plt.ylabel('m')
        locs,labels = plt.xticks()
        # import pdb; pdb.set_trace()
        labels = np.trunc(locs/(24.0*3600.0))
        # locs = np.linspace(-12.0,40,52)
        # labels = range(-12,41)
        plt.xticks(locs,labels)
        
        # Add sea level line
        # t = current_data.t
        plt.hold(True)
        plt.plot([0,0],[0,40],'k-')
        plt.hold(False)
        
    # ========================================================================
    #  Water helper functions
    # ========================================================================
    def b(cd):
        return cd.q[:,:,3] - cd.q[:,:,0]
        
    def extract_eta(h,eta,DRY_TOL=10**-3):
        index = np.nonzero((np.abs(h) < DRY_TOL) + (h == np.nan))
        eta[index[0],index[1]] = np.nan
        return eta
    
    def extract_velocity(h,hu,DRY_TOL=10**-8):
        # u = np.ones(hu.shape) * np.nan
        u = np.zeros(hu.shape)
        index = np.nonzero((np.abs(h) > DRY_TOL) * (h != np.nan))
        u[index[0],index[1]] = hu[index[0],index[1]] / h[index[0],index[1]]
        return u
    
    def eta(cd):
        return extract_eta(cd.q[:,:,0],cd.q[:,:,3])
        
    def water_u(cd):
        # index = np.nonzero(current_data.q[:,:,0] > 1e-6)
        # u = np.zeros(current_data.q[:,:,1].shape)
        # u[index] = current_data.q[index,1] / current_data.q[index,0]
        # return u
        return extract_velocity(cd.q[:,:,0],cd.q[:,:,1])
        # return np.where(abs(current_data.q[:,:,0]) > 10**-16,
        #     current_data.q[:,:,1] / current_data.q[:,:,0],
        #     0.0)
        
    def water_v(cd):
        # index = np.nonzero(current_data.q[:,:,0] > 1e-6)
        # v = np.zeros(current_data.q[:,:,2].shape)
        # v[index] = current_data.q[index,2] / current_data.q[index,0]
        # return u
        return extract_velocity(cd.q[:,:,0],cd.q[:,:,2])
        # return np.where(abs(current_data.q[:,:,0]) > 10**-16,
        #     current_data.q[:,:,2] / current_data.q[:,:,0],
        #     0.0)
        
    def water_speed(current_data):
        u = water_u(current_data)
        v = water_v(current_data)
            
        return np.sqrt(u**2+v**2)
        
    def water_quiver(current_data):
        u = water_u(current_data)
        v = water_v(current_data)
            
        plt.hold(True)
        Q = plt.quiver(current_data.x[::2,::2],current_data.y[::2,::2],
                        u[::2,::2],v[::2,::2])
        max_speed = np.max(np.sqrt(u**2+v**2))
        label = r"%s m/s" % str(np.ceil(0.5*max_speed))
        plt.quiverkey(Q,0.15,0.95,0.5*max_speed,label,labelpos='W')
        plt.hold(False)
            
    def wind_x(cd):
        return cd.q[:,:,4]
    def wind_y(cd):
        return cd.q[:,:,5]
    def wind_speed(cd):
        return np.sqrt(wind_x(cd)**2 + wind_y(cd)**2)

    # ========================================================================
    #  Profile functions
    # ========================================================================
    class PlotProfile(object):
    
        def __init__(self,slice_value = 0.0):
            self.slice_value = slice_value
    
        def slice_index(self,cd):
            if cd.grid.y.lower < self.slice_value < cd.grid.y.upper:
                return int((self.slice_value - cd.grid.y.lower) / cd.dy - 0.5)
            else:
                return None
    
        def bathy_profile(self,current_data):
            index = self.slice_index(current_data)
            if index:
                return current_data.x[:,index], b(current_data)[:,index]
            else:
                return None, None
        
        def surface_profile(self,current_data):
            index = self.slice_index(current_data)
            if index:
                return current_data.x[:,index], eta(current_data)[:,index]
            else:
                return None, None

    # ========================================================================
    #  Plot items
    # ========================================================================
    def add_surface_elevation(plotaxes,bounds=None,plot_type='pcolor'):
        if plot_type == 'pcolor' or plot_type == 'imshow':            
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            # plotitem.plotvar = eta
            plotitem.plot_var = geoplot.surface
            plotitem.imshow_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
            if bounds is not None:
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            plotitem.add_colorbar = True
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1,1,1]
        elif plot_type == 'contour':            
            plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
            plotitem.plot_var = geoplot.surface
            if bounds is None:
                plotitem.contour_levels = [-2.5,-1.5,-0.5,0.5,1.5,2.5]
            # plotitem.contour_nlevels = 21
            # plotitem.contour_min = -2.0
            # plotitem.contour_max = 2.0
            # plotitem.kwargs = {''}
            plotitem.amr_contour_show = [1,1,1]
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1,1,1]
            plotitem.amr_contour_colors = 'k'
            # plotitem.amr_contour_colors = ['r','k','b']  # color on each level
            # plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    
    def add_speed(plotaxes,bounds=None,plot_type='pcolor'):
        if plot_type == 'pcolor' or plot_type == 'imshow':
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            plotitem.plot_var = water_speed
            # plotitem.plot_var = 1
            plotitem.imshow_cmap = plt.get_cmap('PuBu')
            if bounds is not None:
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            plotitem.add_colorbar = True
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1]
        elif plot_type == 'quiver':
            plotitem = plotaxes.new_plotitem(plot_type='2d_quiver')
            plotitem.quiver_var_x = water_u
            plotitem.quiver_var_y = water_v
            plotitem.amr_quiver_show = [4,10,10]
            plotitem.amr_show_key = [True,True,False]
            plotitem.key_units = 'm/s'
            
        elif plot_type == 'contour':
            plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
            plotitem.plot_var = water_speed
            plotitem.kwargs = {'linewidths':1}
            # plotitem.contour_levels = [1.0,2.0,3.0,4.0,5.0,6.0]
            plotitem.contour_levels = [0.5,1.5,3,4.5,6.0]
            plotitem.amr_contour_show = [1,1,1]
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1,1,1]
            plotitem.amr_contour_colors = 'k'
            # plotitem.amr_contour_colors = ['r','k','b']  # color on each level
            # plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
            

    def add_x_velocity(plotaxes,plot_type='pcolor',bounds=None):
        if plot_type == 'pcolor' or plot_type == 'imshow':
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            plotitem.plot_var = water_u
            if bounds is not None:
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            plotitem.add_colorbar = True
            plotitem.imshow_cmap = plt.get_cmap('PiYG')
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1]
        elif plot_type == 'contour':
            pass
    
    def add_y_velocity(plotaxes,plot_type='pcolor',bounds=None):
        if plot_type == 'pcolor' or plot_type == 'imshow':
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            plotitem.plot_var = water_v
            if bounds is not None:
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            plotitem.imshow_cmap = plt.get_cmap('PiYG')
            plotitem.add_colorbar = True
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1]
        elif plot_type == 'contour':
            pass
            
    def add_wind(plotaxes,bounds=None,plot_type='pcolor'):
        if plot_type == 'pcolor' or plot_type == 'imshow':
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            plotitem.plot_var = wind_speed
            plotitem.imshow_cmap = plt.get_cmap('PuBu')
            if bounds is not None:
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            plotitem.add_colorbar = True
            plotitem.amr_imshow_show = [1,1,1]
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1,1,1]
        elif plot_type == 'contour':
            plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
            plotitem.plot_var = wind_speed
            plotitem.contour_nlevels = hurricane_data.max_wind_nest
            plotitem.countour_min = hurricane_data.wind_refine[0]
            plotitem.gridedges_show = 1
        elif plot_type == 'quiver':
            plotitem = plotaxes.new_plotitem(plot_type='2d_quiver')
            plotitem.quiver_var_x = wind_x
            plotitem.quiver_var_y = wind_y
            plotitem.amr_quiver_show = [0,0,1]
            plotitem.amr_quiver_key_show = [True,False,False]
            plotitem.amr_quiver_key_units = 'm/s'
            
    def add_vorticity(plotaxes,bounds=None,plot_type="pcolor"):
        if plot_type == 'pcolor' or plot_type == 'imshow':            
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            plotitem.plot_var = 9
            plotitem.imshow_cmap = plt.get_cmap('PRGn')
            if bounds is not None:
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            plotitem.add_colorbar = True
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1]
            
    def add_land(plotaxes,plot_type='pcolor'):
        if plot_type == 'pcolor':
            plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
            plotitem.show = True
            plotitem.plot_var = geoplot.land
            plotitem.pcolor_cmap = geoplot.land_colors
            plotitem.pcolor_cmin = 0.0
            plotitem.pcolor_cmax = 80.0
            plotitem.add_colorbar = False
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1,1,1]
        elif plot_type == 'contour':            
            plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
            plotitem.plot_var = geoplot.land
            plotitem.contour_nlevels = 40
            plotitem.contour_min = 0.0
            plotitem.contour_max = 100.0
            plotitem.amr_contour_colors = ['g']  # color on each level
            plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
            plotitem.gridlines_show = 0
            plotitem.gridedges_show = 0

    # Limits
    xlimits = [amrdata.xlower,amrdata.xupper]
    ylimits = [amrdata.ylower,amrdata.yupper]
    multilayer_data.eta = eta
    # surface_limits = [-0.15,0.15]
    # speed_limits = [0.0,0.1]
    surface_limits = None
    speed_limits = None
    
    vorticity_limits = [-1.e-2,1.e-2]
    
    # ========================================================================
    #  Surface Elevation
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    
    add_surface_elevation(plotaxes,bounds=surface_limits)
    add_land(plotaxes)
    
    # ========================================================================
    #  Water Speed
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='speed', figno=100)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes

    # Speed
    add_speed(plotaxes,bounds=speed_limits)

    # Land
    add_land(plotaxes)
    
    # X-Velocity
    plotfigure = plotdata.new_plotfigure(name='velocity_x',figno=101)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize':(14,4)}

    # X Velocity
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(121)'
    plotaxes.title = 'X-Velocity'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes

    add_x_velocity(plotaxes,bounds=speed_limits)
    add_land(plotaxes)
    
    # Y Velocity
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(122)'
    plotaxes.title = 'Y-Velocity'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes

    add_y_velocity(plotaxes,bounds=speed_limits)
    add_land(plotaxes)
    
    
    # ========================================================================
    #  Wind speed
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='wind',figno=2)
    plotfigure.show = True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Wind"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    
    # Wind
    add_wind(plotaxes)
    add_land(plotaxes)

    # ========================================================================
    #  Profile Plots
    # ========================================================================
    # Profile variables
            
    def profile_afteraxes(current_data):
        day_figure_title(current_data)
        loc,label = plt.xticks()
        label = loc/1.e3
        plt.xticks(loc,label)
        plt.xlabel('km')
        if current_data.plotaxes.title == 'Wind':
            plt.ylabel('m/s')
        else:
            plt.ylabel('m')
            
        t = current_data.t
        # Hurricane eye
        x = t * hurricane_data.hurricane_velocity[0] + hurricane_data.R_eye_init[0]
        plt.hold(True)
        plt.plot(x,0.0,'r+')
        plt.hold(False)
            
    def remove_labels_profile(cd,direction='x'):
        plt.hold(True)
        if direction == 'x':
            plt.xlabel('')
            locs,labels = plt.xticks()
            # labels = np.flipud(locs)/1.e3
            labels = ['' for i in xrange(len(locs))]
            plt.xticks(locs,labels)
            plt.ylabel('m')
        elif direction == 'y':
            plt.ylabel('')
            locs,labels = plt.yticks()
            # labels = np.flipud(locs)/1.e3
            labels = ['' for i in xrange(len(locs))]
            plt.yticks(locs,labels)
            plt.xlabel('m')
        plt.hold(False)
        
    def labels_profile(cd,direction='x'):
        if direction == 'x':
            loc,label = plt.xticks()
            label = loc/1.e3
            plt.xticks(loc,label)
            plt.xlabel('km')
            if cd.plotaxes.title == 'Wind':
                plt.ylabel('m/s')
            else:
                plt.ylabel('m')
        elif direction == 'y':
            loc,label = plt.yticks()
            label = loc/1.e3
            plt.yticks(loc,label)
            plt.ylabel('km')
            if cd.plotaxes.title == 'Wind':
                plt.xlabel('m/s')
            else:
                plt.xlabel('m')
        
    def profile_afteraxes(current_data):
        day_figure_title(current_data)
        labels_profile(current_data)
        # bathy_ref_lines_profile(current_data,surface_limits)
    
    
    plotfigure = plotdata.new_plotfigure(name='profile', figno=4)
    plotfigure.show = False
        
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Profiles'
    plotaxes.xlimits = xlimits
    # plotaxes.ylimits = surface_limits
    plotaxes.afteraxes = profile_afteraxes
    
    profile_plot = PlotProfile(0.0)
    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
    plotitem.map_2d_to_1d = profile_plot.surface_profile
    plotitem.amr_plotstyle = ['-','-.','+','x','.']
    plotitem.color = 'b'#(0.2,0.8,1.0)
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = profile_plot.bathy_profile
    plotitem.amr_plotstyle = ['-','-.','+','x','.']  
    plotitem.color = 'k'
        
    # ========================================================================
    #  Bathy Profile
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='bathy_profile',figno=20)
    plotfigure.show = False
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.title = "Bathymetry Profile"
    plotaxes.scaled = 'equal'
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = b
    plotitem.imshow_cmap = plt.get_cmap('earth')
    plotitem.imshow_cmin = -3300
    plotitem.imshow_cmax = 100.0
    plotitem.add_colorbar = True
    plotitem.amr_imshow_show = [1,1,1]
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1,1,1]
    plotitem.show = True
    
    # ========================================================================
    # Figure for grids alone
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='grids', figno=11)
    plotfigure.show = False
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.xlimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.title = 'grids'
    plotaxes.afteraxes = pcolor_afteraxes
    plotaxes.scaled = True
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_grid')
    # plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_grid_bgcolor = ['blue','red','green','cyan','yellow']
    plotitem.amr_gridlines_show = [1,1,0,0,0,0]   
    plotitem.amr_gridedges_show = 1
    
    # ========================================================================
    # Figures for momentum
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='x-momentum', figno=13)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'X-Velocity'
    plotaxes.scaled = True
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.afteraxes = pcolor_afteraxes
    
    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.surface
    plotitem.plot_var = water_u
    plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    # plotitem.pcolor_cmin = -1.e-10
    # plotitem.pcolor_cmax = 1.e-10
    # plotitem.pcolor_cmin = -2.5 # -3.0
    # plotitem.pcolor_cmax = 2.5 # 3.0
    plotitem.add_colorbar = True
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1,1,1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.show = True
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 80.0
    plotitem.add_colorbar = False
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1,1,1]
    
    plotfigure = plotdata.new_plotfigure(name='y-momentum', figno=14)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Y-Velocity'
    plotaxes.scaled = True
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.afteraxes = pcolor_afteraxes
    
    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.surface
    plotitem.plot_var = water_v
    plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    # plotitem.pcolor_cmin = -1.e-10
    # plotitem.pcolor_cmax = 1.e-10
    # plotitem.pcolor_cmin = -2.5 # -3.0
    # plotitem.pcolor_cmax = 2.5 # 3.0
    plotitem.add_colorbar = True
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1,1,1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.show = True
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 80.0
    plotitem.add_colorbar = False
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1,1,1]
    
    # ========================================================================
    #  Contour plot for surface
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='contour_surface',figno=15)
    plotfigure.show = False
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.afteraxes = contour_afteraxes
    
    # Surface
    add_surface_elevation(plotaxes,plot_type='contour')
    
    # Land
    add_land(plotaxes,plot_type='contour')
    
    # ========================================================================
    #  Contour plot for speed
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='contour_speed',figno=16)
    plotfigure.show = False
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Current'
    plotaxes.scaled = True
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.afteraxes = contour_afteraxes
    
    # Surface
    add_surface_elevation(plotaxes,plot_type="contour")
    
    # Land
    add_land(plotaxes,plot_type='contour')
    
    # ========================================================================
    #  Vorticity Plot
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='vorticity',figno=17)
    plotfigure.show = False
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Vorticity"
    plotaxes.scaled = True
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.afteraxes = pcolor_afteraxes
    
    # Vorticity
    add_vorticity(plotaxes)

    # Land
    add_land(plotaxes)
    
    
    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')
    plotfigure.show = False
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.0,40.0*3600.0]
    # plotaxes.ylimits = [0,150.0]
    plotaxes.ylimits = [-3.0, 3.0]
    plotaxes.title = 'Surface'
    plotaxes.afteraxes = gauge_afteraxes

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'r-'

    # # Plot topo as green curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = gaugetopo
    # plotitem.plotstyle = 'g+'

    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                    # create html files of plots?
    plotdata.latex = False                   # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
