#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 12:17:12 2020

@author: Mark Richmond
"""
from scipy.interpolate import CubicSpline
from pylab import *
import cmath
from scipy import optimize

'''
This program optimises the anti-reflection refractive index and layer
thickness for a given normal incident wavelength. It is set up for an air
medium and a BK7 substrate.
'''

#Fit cubic spline to real and imaginary parts of refractive index data
w0, n0, k0 = loadtxt("Air.txt", skiprows=1, unpack=True)
N0 = CubicSpline(w0, n0)
K0 = CubicSpline(w0, k0)

w2, n2, k2 = loadtxt("BK7.txt", skiprows=1, unpack=True)
N2 = CubicSpline(w2, n2)
K2 = CubicSpline(w2, k2)

#This is the wavelength in nm of incident light
wavelength = 500

'''
Function (calculating reflectance) to be minimised through the optimisation
of 2 variables: coating refractive index and coating thickness.
'''

def z(params):
    N1, layer_thickness = params
    Refractive_real = array([N0(wavelength),N1,N2(wavelength)])
    Refractive_imaginary = array([K0(wavelength),0,K2(wavelength)])*(0+1j)
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

#Initial guess as starting point for optimization.
initial_guess = [1.2, 100]
result = optimize.minimize(z, initial_guess)
if result.success:
    fitted_params = result.x
    print ('Reflectance minimsed at [Refractive index, Layer thickness (nm)]')
    print (fitted_params)
    print ('Refractive index of air')
    print (N0(wavelength))
    print ('Refractive index of BK7')
    print (N2(wavelength))
    print ('Theoretical prediction of optimal refractive index, sqrt(n0n2)')
    print (sqrt(N0(wavelength)*N2(wavelength)))    
    print ('Optical thickness of coating (Number of wavelengths)')
    print ((fitted_params[0]*fitted_params[1])/wavelength)
else:
    raise ValueError(result.message)

