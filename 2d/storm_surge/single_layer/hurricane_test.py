#!/usr/bin/env python
# encoding: utf-8
r"""
Brief description

In depth description

:Authors:
    Kyle Mandli (2009-12-15) Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import sys
import os

import numpy as np
from pylab import *

from pyclaw.solution import Solution

if not os.path.exists('./hurricane'):
    os.mkdir('./hurricane')

for i in xrange(11):
    sol = Solution(i)
    [X,Y] = sol.grids[0].p_center
    clf()
    # subplot(2,1,1)
    quiver(X,Y,sol.q[:,:,0],sol.q[:,:,1])
    axis('equal')
    title('Wind speed, t = %s h' % str(sol.t/(60**2)))
        
    # subplot(2,1,2)    # 
        # pcolor(X,Y,sol.q[:,:,2])
        # axis('equal')
        # title('Pressure Field, t = %s h' % str(sol.t/(60**2)))
        # colorbar()
    
    print "Saving frame ", i
    savefig('hurricane/frame%s.png' % str(i).zfill(4))
    
