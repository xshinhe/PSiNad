# distutils: language = c++
# cython: language_level = 3

"""
test.pyx: Utility functions that are compiled with Cython for speed

This is a demo file.
Authors: Xin He
Contributors:
"""

from __future__ import absolute_import
__author__ = "Xin He"
__version__ = "1.0"

# from heapq import heappush, heappop

def fun1_add(a, b):
    return a + b
