
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

    amrdata = Data(os.path.join(plotdata.rundir,'amr2ez.data'))
    hurricane_data = Data(os.path.join(plotdata.rundir,'hurricane.data'))
    plotdata.clearfigures()
    plotdata.clear_frames = False
    plotdata.clear_figs = True
    
    plotdata.save_frames = False
    
    # ========================================================================
    #  Generic helper functions
    # ========================================================================
    def pcolor_afteraxes(current_data):
        eye_location(current_data)
        hour_figure_title(current_data)
        m_to_km_labels()
        # wind_contours(current_data)
        bathy_ref_lines(current_data)
        
    def contour_afteraxes(current_data):
        eye_location(current_data)
        hour_figure_title(current_data)
        # gauge_locations(current_data)
        m_to_km_labels()
        # plt.hold(True)
        # pos = -80.0 * (23e3 / 180) + 500e3 - 5e3
        # plt.plot([pos,pos],[-300e3,300e3],'b',[pos-5e3,pos-5e3],[-300e3,300e3],'y')
        # plt.hold(False)
        # wind_contours(current_data)
        # bathy_ref_lines(current_data)
        
    def hour_figure_title(current_data):
        t = current_data.t
        title = current_data.plotaxes.title
        plt.title('%s at time t = %s h' % (title,str(t/3600.0)))

    def m_to_km_labels(current_data=None):
        plt.xlabel('km')
        plt.ylabel('km')
        locs,labels = plt.xticks()
        labels = locs/1.e3
        plt.xticks(locs,labels)
        locs,labels = plt.yticks()
        labels = locs/1.e3
        plt.yticks(locs,labels)

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
        plt.xlabel('t (hours)')
        plt.ylabel('m')
        locs,labels = plt.xticks()
        # import pdb; pdb.set_trace()
        labels = np.trunc(locs/3600.0)
        # locs = np.linspace(-12.0,40,52)
        # labels = range(-12,41)
        plt.xticks(locs,labels)
        
        # Add sea level line
        # t = current_data.t
        plt.hold(True)
        plt.plot([0,0],[0,40],'k-')
        plt.hold(False)

        
    # ========================================================================
    #  Hurricane related helper functions
    # ========================================================================
    # Hurricane eye location
    def eye_location(current_data):
        t = current_data.t
        # Hurricane eye
        x = t * hurricane_data.hurricane_velocity[0] + hurricane_data.R_eye_init[0]
        y = t * hurricane_data.hurricane_velocity[1] + hurricane_data.R_eye_init[1]
        
        plt.hold(True)
        plt.plot(x,y,'rD')
        plt.hold(False)
        
    def hurricane_afteraxes(current_data):
        eye_location(current_data)
        hour_figure_title(current_data)
        m_to_km_labels()
        
    def bathy_ref_lines(current_data):
        plt.hold(True)
        y = [amrdata.ylower,amrdata.yupper]
        # for ref_line in [352e3,410e3,477e3]:
        for ref_line in [450e3]:
            plt.plot([ref_line,ref_line],y,'y--')
        plt.hold(False)
        
    def hurricane_wind(current_data):
        if current_data.level == 1:
            t = current_data.t
            u = current_data.q[:,:,8]
            v = current_data.q[:,:,9]
            plt.hold(True)
            Q = plt.quiver(current_data.x[::3,::3],current_data.y[::3,::3],
                        u[::3,::3],v[::3,::3])
            # plt.quiverkey(Q,0.5,0.5,50,r'$50 \frac{m}{s}$',labelpos='W',
            #                 fontproperties={'weight':'bold'})
            plt.hold(False)
            
    def wind_speed(current_data):
        return np.sqrt(current_data.q[:,:,4]**2 + current_data.q[:,:,5]**2)
        
    def wind_contours(current_data):
        plt.hold(True)
        w = wind_speed(current_data)
        max_w = np.max(np.max(w))
        levels = [0.0,0.25*max_w,0.5*max_w,0.75*max_w,max_w*0.999]
        C = plt.contour(current_data.x,current_data.y,w,levels)
        plt.clabel(C,inline=1)
        plt.hold(False)
        
    # ========================================================================
    #  Water helper functions
    # ========================================================================
    def b(cd):
        return cd.q[:,:,3] - cd.q[:,:,0]
        
    def eta(cd):
        return cd.q[:,:,3]
    
    def water_u(current_data):
        # index = np.nonzero(current_data.q[:,:,0] > 1e-6)
        # u = np.zeros(current_data.q[:,:,1].shape)
        # u[index] = current_data.q[index,1] / current_data.q[index,0]
        # return u
        return np.where(abs(current_data.q[:,:,0]) > 10**-16,
            current_data.q[:,:,1] / current_data.q[:,:,0],
            0.0)
        
    def water_v(current_data):
        # index = np.nonzero(current_data.q[:,:,0] > 1e-6)
        # v = np.zeros(current_data.q[:,:,2].shape)
        # v[index] = current_data.q[index,2] / current_data.q[index,0]
        # return u
        return np.where(abs(current_data.q[:,:,0]) > 10**-16,
            current_data.q[:,:,2] / current_data.q[:,:,0],
            0.0)
        
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
    
    # ========================================================================
    #  Surface Elevation
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.afteraxes = pcolor_afteraxes
    
    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = geoplot.surface
    # plotitem.plot_var = 0
    # plotitem.imshow_cmin = -2.5
    # plotitem.imshow_cmax = 2.5
    # plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.imshow_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    # plotitem.pcolor_cmin = -1.e-2
    # plotitem.pcolor_cmax = 1.e-2
    # plotitem.pcolor_cmin = -1.0
    # plotitem.pcolor_cmax = 1.0
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
    #  Water Speed
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='speed', figno=1)
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents'
    plotaxes.scaled = True
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    
    def pcolor_afteraxes(current_data):
        eye_location(current_data)
        hour_figure_title(current_data)
        # bathy_ref_lines(current_data)
        m_to_km_labels()
    plotaxes.afteraxes = pcolor_afteraxes

    # Speed
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = water_speed
    # plotitem.plot_var = 1
    plotitem.imshow_cmap = plt.get_cmap('PuBu')
    # plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    # plotitem.pcolor_cmin = 0.0
    # plotitem.pcolor_cmax = 1e-1 # 6.0
    # plotitem.pcolor_cmin = 0.0
    # plotitem.pcolor_cmax = 5.0 # 6.0
    # plotitem.imshow_cmin = 0.0
    # plotitem.imshow_cmax = 0.0015
    plotitem.add_colorbar = True
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1]
    # plotitem.aftergrid = water_quiver
    
    # Velocity vectors
    plotitem = plotaxes.new_plotitem(plot_type='2d_quiver')
    plotitem.quiver_var_x = water_u
    plotitem.quiver_var_y = water_v
    plotitem.amr_quiver_show = [4,10,10]
    plotitem.amr_show_key = [True,True,False]
    plotitem.key_units = 'm/s'
    plotitem.show = False

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 80.0
    plotitem.add_colorbar = False
    plotitem.amr_gridlines_show = [0,0,0]
    
    #-----------------------------------------
    # Hurricane forcing
    #-----------------------------------------
    # Pressure field
    plotfigure = plotdata.new_plotfigure(name='pressure', figno=2)
    plotfigure.show = hurricane_data.pressure_src and True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.title = "Pressure Field"
    plotaxes.afteraxes = hurricane_afteraxes
    plotaxes.scaled = True
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    def pressure_var(current_data):
        return current_data.q[:,:,6] / 100.0
    plotitem.plot_var = pressure_var
    plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    plotitem.pcolor_cmin = 954
    plotitem.pcolor_cmax = 1002
    plotitem.add_colorbar = True
    plotitem.gridlines_show = 0
    plotitem.gridedges_show = 1
    plotitem.amr_pcolor_show = [1,0,0]
    plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 80.0
    plotitem.add_colorbar = False
    plotitem.amr_gridlines_show = [0,0,0]

    # Pressure gradient plots - x
    plotfigure = plotdata.new_plotfigure(name='pressure_x', figno=30)
    plotfigure.show = hurricane_data.pressure_src and True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.title = "Pressure Gradient Field - x"
    plotaxes.afteraxes = hurricane_afteraxes
    plotaxes.scaled = True
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    def pressure_var(current_data):
        return current_data.q[:,:,7] / 100.0
    plotitem.plot_var = pressure_var
    plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    # plotitem.pcolor_cmin = 954
    # plotitem.pcolor_cmax = 1002
    plotitem.add_colorbar = True
    plotitem.gridlines_show = 0
    plotitem.gridedges_show = 1
    plotitem.amr_pcolor_show = [1,0,0]
    plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    
    # Pressure gradient plots - y
    plotfigure = plotdata.new_plotfigure(name='pressure_y', figno=31)
    plotfigure.show = hurricane_data.pressure_src and True
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.title = "Pressure Gradient Field - y"
    plotaxes.afteraxes = hurricane_afteraxes
    plotaxes.scaled = True
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    def pressure_var(current_data):
        return current_data.q[:,:,8] / 100.0
    plotitem.plot_var = pressure_var
    plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    # plotitem.pcolor_cmin = 954
    # plotitem.pcolor_cmax = 1002
    plotitem.add_colorbar = True
    plotitem.gridlines_show = 0
    plotitem.gridedges_show = 1
    plotitem.amr_pcolor_show = [1,0,0]
    plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    
    
    # ========================================================================
    # Wind field
    plotfigure = plotdata.new_plotfigure(name='wind',figno=3)
    plotfigure.show = hurricane_data.wind_src and False
    
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.ylimits = [amrdata.ylower,amrdata.yupper]
    plotaxes.title = "Wind Field"
    plotaxes.afteraxes = hurricane_afteraxes
    plotaxes.scaled = True
    
    # Quiver
    plotitem = plotaxes.new_plotitem(plot_type='2d_quiver')
    plotitem.quiver_var_x = 4
    plotitem.quiver_var_y = 5
    plotitem.amr_quiver_show = [0,0,1]
    plotitem.amr_quiver_key_show = [True,False,False]
    plotitem.amr_quiver_key_units = 'm/s'
    plotitem.show = False
    
    # Pcolor
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = wind_speed
    plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    plotitem.pcolor_cmin = 0
    plotitem.pcolor_cmax = 55
    plotitem.add_colorbar = True
    plotitem.amr_pcolor_show = [1,1,1]
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1,1,1]
    plotitem.show = True
    
    # Contour
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = wind_speed
    plotitem.contour_nlevels = hurricane_data.max_wind_nest
    plotitem.countour_min = hurricane_data.wind_refine[0]
    plotitem.gridedges_show = 1
    plotitem.show = False 
    
    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 80.0
    plotitem.add_colorbar = False
    plotitem.amr_gridlines_show = [0,0,0]

    # ========================================================================
    #  Profile Plots
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='profile', figno=4)
    plotfigure.show = False
    
    def after_axes_profile(current_data):
        hour_figure_title(current_data)
        loc,lable = plt.xticks()
        lable = loc/1.e3
        plt.xticks(loc,lable)
        plt.xlabel('km')
        if current_data.plotaxes.title == 'Wind':
            plt.ylabel('m/s')
        else:
            plt.ylabel('m')
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.title = "Wind"
    plotaxes.afteraxes = after_axes_profile
    
    # Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = geoplot.surface
    def sea_surface_profile(current_data):
    #     if current_data.level == 1:
    #         q = current_data.var
    #         return current_data.x[:,0],q[:,0]
        return None
    plotitem.map_2d_to_1d = sea_surface_profile
    plotitem.plotstyle = '-k'
    plotitem.show = False
    
    # Speed
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = water_speed
    def sea_speed_profile(current_data):
        if current_data.level == 1:
            q = current_data.var
            x = current_data.x
            return x[:,0],q[:,0]
        return None
    plotitem.map_2d_to_1d = sea_speed_profile
    plotitem.plotstyle = '-k'
    plotitem.show = False
    
    # Wind
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    def wind_speed(current_data):
        return np.sqrt(current_data.q[:,:,4]**2 + current_data.q[:,:,5]**2)
    def wind_speed_profile(current_data):
        q = current_data.var
        index = np.floor(amrdata.my / 2)
        return current_data.x[:,index],q[:,index]
    plotitem.plot_var = wind_speed
    plotitem.map_2d_to_1d = wind_speed_profile
    plotitem.plotstyle = '-k'
    plotitem.show = False
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = [amrdata.xlower,amrdata.xupper]
    plotaxes.title = "Bathymetery"
    plotaxes.afteraxes = after_axes_profile
    plotitem.amr_color=['r','b','g']
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')    
    def bathy(current_data):
        return current_data.q[:,:,3] - current_data.q[:,:,0]
    plotitem.plot_var = bathy
    def bath_profile(current_data):
        q = current_data.var
        index = np.floor(amrdata.my / 2)
        return current_data.x[:,index],q[:,index]
    plotitem.map_2d_to_1d = bath_profile
    plotitem.plotstyle = '-k'
    plotitem.show = False
        
        
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
    #  Scatter plot of surface for radially symmetric
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Scatter', figno=12)
    plotfigure.show = False
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.,400e3]
    # plotaxes.ylimits = [-2.5,]
    plotaxes.title = 'Scatter plot of surface'
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = geoplot.surface
    def q_vs_radius(current_data):
        q = current_data.var
        return current_data.x[:,20],q[:,20]
    plotitem.map_2d_to_1d = q_vs_radius
    plotitem.plotstyle = 'o-'
    plotitem.amr_color=['b','r','g']
    
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = geoplot.land
    def q_vs_radius(current_data):
        q = current_data.var
        return current_data.x[:,20],q[:,20]
    plotitem.map_2d_to_1d = q_vs_radius
    plotitem.plotstyle = 'o-'
    plotitem.amr_color=['r','r','g']
    plotaxes.afteraxes = "pylab.legend(['Level 1','Level 2'])"
    
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
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.surface
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
    plotitem.show = True 
    
    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.land
    plotitem.contour_nlevels = 40
    plotitem.contour_min = 0.0
    plotitem.contour_max = 100.0
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.gridlines_show = 0
    plotitem.gridedges_show = 0
    plotitem.show = True
    
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
    plotitem.show = True 
    
    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.land
    plotitem.contour_nlevels = 40
    plotitem.contour_min = 0.0
    plotitem.contour_max = 100.0
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.gridlines_show = 0
    plotitem.gridedges_show = 0
    plotitem.show = True
    
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
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = 9
    plotitem.imshow_cmap = plt.get_cmap('PRGn')
    # plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    # plotitem.pcolor_cmin = 0.0
    # plotitem.pcolor_cmax = 6.0
    plotitem.imshow_cmin = -1.e-2
    plotitem.imshow_cmax = 1.e-2
    plotitem.add_colorbar = True
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 80.0
    plotitem.add_colorbar = False
    plotitem.amr_gridlines_show = [0,0,0]
    
    
    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')
    plotfigure.show = True
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

    
