#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing test cases

:Authors:
    Kyle Mandli (2011-05-07) Initial version
"""
# ============================================================================
#      Copyright (C) 2011 Kyle Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import sys
import numpy as np

# Tests
tests = [{"name":"perpendicular","velocity":5.0, "angle": 0.00 * np.pi, "eye":(0.0,0.0)},
         {"name":"p45angle","velocity":5.0, "angle": 0.25 * np.pi, "eye":(200e3,100e3)},
         {"name":"n45angle","velocity":5.0, "angle": -0.25 * np.pi, "eye":(200e3,-100e3)},
         {"name":"parallel_down","velocity":5.0, "angle": -0.50 * np.pi, "eye":(425e3,100e3)},
         {"name":"parallel_up","velocity":5.0, "angle": 0.50 * np.pi, "eye":(425e3,-100e3)},
         {"name":"static","velocity":5.0,"angle":0.0, "eye":(0.0,0.0)}]

def print_tests():
    for (i,test) in enumerate(tests):
        v =(test["velocity"] * np.cos(test["angle"]),test["velocity"] * np.sin(test["angle"]))
        print "Test %s: %s" % (i,test['name'])
        print "   velocity = (%s,%s) - angle = %s" % (v[0],v[1],test['angle'])
        print "   eye = (%s,%s)" % (test['eye'][0],test['eye'][1])

if __name__ == "__main__":
    sys.exit(print_tests())
