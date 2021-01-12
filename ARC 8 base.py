#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 11:32:34 2020

@author: Mark Richmond
"""
from scipy.interpolate import CubicSpline
from pylab import *
import cmath

'''
This program plots the reflectance spectrum for a single layer film on a
substrate in a given medium at normal incidence. It is set up for a MgF2 layer
on a substrate of BK7 in air.
'''

#Fit cubic spline to real and imaginary parts of refractive index data
w0, n0, k0 = loadtxt("Air.txt", skiprows=1, unpack=True)
N0 = CubicSpline(w0, n0)
K0 = CubicSpline(w0, k0)

w1, n1, k1 = loadtxt("MgF2.txt", skiprows=1, unpack=True)
N1 = CubicSpline(w1, n1)
K1 = CubicSpline(w1, k1)

w2, n2, k2 = loadtxt("BK7.txt", skiprows=1, unpack=True)
N2 = CubicSpline(w2, n2)
K2 = CubicSpline(w2, k2)

'''
Function to calculate the overall transfer matrix, Z, and to return the
reflectance as calculated from the appropriate matrix elements. Equations for
s-polarization are used but results will be the same for p-polarization at
normal incidence. Single layer between medium and substrate.
'''

def z(wavelength,layer_thickness):
    Refractive_real = array([N0(wavelength),N1(wavelength),N2(wavelength)])
    Refractive_imaginary = array([K0(wavelength),K1(wavelength),K2(wavelength)])*(0+1j)
    R = Refractive_real+Refractive_imaginary
    K=((2*pi)/wavelength)*R
    T=zeros((2,2),complex)
    T[0,0]= (0.5)*(sqrt(K[1]/K[0])+sqrt(K[0]/K[1]))
    T[0,1]= (0.5)*(sqrt(K[1]/K[0])-sqrt(K[0]/K[1]))
    T[1,0]=T[0,1]
    T[1,1]=T[0,0]
    P=zeros((2,2),complex)
    P[0,0]=cos(K[1]*layer_thickness)+((sin(K[1]*layer_thickness))*(0+1j))
    P[1,1]=cos(K[1]*layer_thickness)-((sin(K[1]*layer_thickness))*(0+1j))
    X=matmul(P,T)
    Tn=zeros((2,2),complex)
    Tn[0,0]= (0.5)*(sqrt(K[2]/K[1])+sqrt(K[1]/K[2]))
    Tn[0,1]= (0.5)*(sqrt(K[2]/K[1])-sqrt(K[1]/K[2]))
    Tn[1,0]=Tn[0,1]
    Tn[1,1]=Tn[0,0]
    Z=matmul(Tn,X)
    r=-(Z[1,0]/Z[1,1])
    return abs(r)**2
'''
plots the reflectance for wavelengths between 350 and 1000nm. At the moment
it is setup for 90nm layer thickness. To plot a more detailed graph reduce the
number taken as the third argument in the arange command below.
'''

axis([350,1000,0.01,0.03])
for i in arange(350,1000,0.1):
    
    plt.plot(i,z(i,90), 'b.',ms=1)
xlabel('Wavelength (nm)')
ylabel('Reflectance')
grid()