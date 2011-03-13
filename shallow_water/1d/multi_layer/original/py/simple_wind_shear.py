#!/usr/bin/env python
# encoding: utf-8
r"""
Simple wind shear example utilizing the multi-layer shallow water solver

Based off of the solver in:

E. Audusse, M. Bristeau, B. Perthame, and J. Sainte-Marie. A multilayer 
saint-venant system with mass exchanges for shallow water ﬂows. Mathematical 
Modelling and Numerical Analysis, 2009.

:Authors:
    Kyle T. Mandli (2009-9-2):  Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import numpy as np

from pyclaw.controller import Controller
from pyclaw.solution import Dimension, Grid, Solution
from pyclaw.evolve.clawpack import ClawSolver1D

from rp_multi_layer_swe import rp_multi_layer_swe_1d_linearized

# Parameters
g = 1.0
layers = 2

# Setup controller
claw = Controller()

# Setup initial condition
x = Dimension('x',0.0,1.0,200)
claw.solution = Solution(x)
claw.solution.aux_global['layers'] = layers
claw.solution.aux_global['g'] = g
claw.solution.meqn = layers + 1
claw.solution.q = np.zeros((x.n,claw.solution.meqn))

# Setup solver
claw.solver = ClawSolver1D()
claw.solver.rp = rp_multi_layer_swe_1d_linearized

# Run
claw.run()

# Plot