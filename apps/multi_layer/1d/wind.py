#!/usr/bin/env python

r"""Module containing functions pertaining to wind field calculations"""

# =============================================================================
#  Wind functions
def set_constant_wind(state,w):
    r"""Set a contant wind field"""
    state.aux[1,:] = w

# No wind field convenience function
set_no_wind = lambda state: set_constant_wind(state,0.0)

def set_oscillatory_wind(state,A,N,omega,t_length):
    r"""Set wind to an oscillating wind field"""
    x = state.grid.x.centers
    t = state.t
    L = state.grid.x.upper - state.grid.x.lower
    
    state.aux[1,:] = A * np.sin(np.pi * N * x / L) \
                       * np.sin(2.0 * np.pi * omega / t_length * t)

def wind_drag(wind_speed):
    r"""Calculate wind drag coefficient"""
    if wind_speed <= 11.0:
        wind_drag = 1.20
    elif wind_speed > 11.0 and wind_speed <= 25.0:
        wind_drag = 0.49 + 0.065 * wind_speed
    else:
        wind_drag = 0.49 + 0.065 * 25.0
        
    return wind_drag * 10.e-3
    