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

import test_suites

tests = []

# Two Layer Tests
tests.append(TwoLayerBaseTest(5.0,0.0,(0.0,0.0),5))
tests.append(TwoLayerBaseTest(5.0,0.25*np.pi,(200e3,100e3),5))
tests.append(TwoLayerBaseTest(5.0,-0.25*np.pi,(200e3,-100e3),5))
tests.append(TwoLayerBaseTest(5.0,-0.50*np.pi,(425e3,100e3),5))
tests.append(TwoLayerBaseTest(5.0, 0.50*np.pi,(425e3,-100e3),5))
tests.append(TwoLayerBaseTest(0.0,0.0,(0.0,0.0),5))

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'all':
            tests_to_be_run = tests
        else:
            tests_to_be_run = []
            for test in sys.argv[1:]:
                tests_to_be_run.append(tests[int(test)])
            
        test_runs.run_tests(tests,parallel=True)

    else:
        test_runs.print_tests(tests)