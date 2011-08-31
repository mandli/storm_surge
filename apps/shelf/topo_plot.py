#!/usr/bin/env python
# encoding: utf-8
r"""
Plot profile and overview plot of topography

:Authors:
    Kyle Mandli (2010-07-09) Initial version
"""

import os

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('./topo.data')
X = data[:,0].reshape((5,5))
Y = data[:,1].reshape((5,5))
Z = data[:,2].reshape((5,5))

print Z

plt.pcolor(Z)
plt.colorbar()
plt.show()