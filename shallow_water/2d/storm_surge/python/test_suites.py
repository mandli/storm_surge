#!/usr/bin/env python
# encoding: utf-8
r"""
Runs tests for angles relative to a continental shelf

:Authors:
    Kyle Mandli (2011-04-20) Initial version
"""
# ============================================================================
#      Copyright (C) 2011 Kyle Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import os
import sys

import numpy as np

import test_runs

class SingleLayerBaseTest(test_runs.Test):
    
    def __init__(self,velocity=5.0,angle=0.0,eye=(0.0,0.0),mxnest=5):
        self.type = "storm_surge"
        self.name = "single_layer"
        self.setplot = "setplot"
        
        self.run_data = setrun.setrun()
        self.ml_data = setrun.set_multilayer_data()
        self.hurricane_data = setrun.set_hurricane_data()
        
        self.run_data.mxnest = -mxnest
        self.hurricane_velocity = (velocity * np.cos(angle),velocity * np.sin(angle))
        self.R_eye_init = eye
        
        self.prefix = "sl_angle%s_m%s" % (int(angle * 180.0 / np.pi),mxnest)

class TwoLayerBaseTest(test_runs.TestML2D,SingleLayerBaseTest):
    
    def __init__(self,velocity=5.0,angle=0.0,eye=(0.0,0.0),mxnest=5):
        super(TwoLayerBaseTest,self).__init__()
        
        self.name = "two_layer"
        
        self.prefix = "ml_angle%s_m%s" % (int(angle * 180.0 / np.pi),mxnest)

