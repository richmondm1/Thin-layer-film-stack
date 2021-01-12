#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 16:17:30 2020

@author: Mark Richmond
"""

from scipy.interpolate import CubicSpline
from scipy.optimize import leastsq
from pylab import *

'''
This program produces a cubic spline for a given material and returns the
interpolated values of the real and imaginary components of the refractive
index at a given wavelength. It also plots the real and imaginary components
of the refractive index. A file of spectral refractive index data points is
required as an input. Ta205 data is processed here but the input file could be
altered for other materials used in this group of programs: MgF2, Air, Au, BK7,
AlAs.
'''

w, n, k = loadtxt("Ta2O5.txt", skiprows=1, unpack=True)
W = linspace(min(w),max(w),50000)
N = CubicSpline(w, n)
K = CubicSpline(w, k)
print (N(633))
print (K(633))


plot(w,n,'b,')
plot(W,N(W),'b-')
plot(w,k,'r,')
plot(W,K(W),'r-')
grid()
xlabel('Wavelength (nm)')
ylabel('Refractive index')

