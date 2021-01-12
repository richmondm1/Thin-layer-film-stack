#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 13:14:37 2021

@author: Mark Richmond
"""
from scipy.interpolate import CubicSpline
from pylab import *
import cmath
from numpy.lib.scimath import arcsin

'''
This program calculates reflectance, transmittance and absorption from a
single layer film on a substrate in a medium. It is setup for a gold film on
BK7 in air. Normal incidence. Currently configured to plot reflectance graph
but transmission or absorption plotting can alternatively be selected by
activating the relevant sections at the bottom of the code.
'''


#Fit cubic spline to real and imaginary parts of refractive index data
w0, n0, k0 = loadtxt("Air.txt", skiprows=1, unpack=True)
N0 = CubicSpline(w0, n0)
K0 = CubicSpline(w0, k0)

w1, n1, k1 = loadtxt("Au.txt", skiprows=1, unpack=True)
N1 = CubicSpline(w1, n1)
K1 = CubicSpline(w1, k1)

w2, n2, k2 = loadtxt("BK7.txt", skiprows=1, unpack=True)
N2 = CubicSpline(w2, n2)
K2 = CubicSpline(w2, k2)

'''
Function to calculate the overall transfer matrix, Z, and to return the
reflectance, transmittance and absorption as calculated from the appropriate
matrix elements. Equations for s-polarization are used. Single layer between
medium and substrate.
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
    t=Z[0,0]+(Z[0,1]*r)
    Reflectance=abs(r)**2
    Transmittance=abs(t)**2
    Absorption=1-Reflectance-Transmittance
    return [Reflectance, Transmittance, Absorption] 

'''
Plots the reflectance for wavelengths between 350 and 1000nm. To plot a more
detailed graph reduce the number taken as the third argument in the arange 
commands below.
'''

axis([350,1000,0,1])
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,10)[0], 'b.',label="10nm" if i==350 else "",ms=1)
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,50)[0], 'r.',label="50nm" if i==350 else "",ms=1)
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,100)[0], 'g.',label="100nm" if i==350 else "",ms=1)
xlabel('Wavelength (nm)')
ylabel('Reflectance')
plt.legend(loc="upper left",markerscale=4)
grid()

'''
#plots the transmittance for wavelengths between 350 and 1000nm.
#please comment out the plots above and below if you wish to see this plot.
axis([350,1000,0,1])
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,10)[1], 'b.',label="10nm" if i==350 else "",ms=1)
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,50)[1], 'r.',label="50nm" if i==350 else "",ms=1)
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,100)[1], 'g.',label="100nm" if i==350 else "",ms=1)
xlabel('Wavelength (nm)')
ylabel('Transmittance')
plt.legend(loc="upper left",markerscale=4)
grid()

#plots the absorption for wavelengths between 350 and 1000nm.
#please comment out the 2 plots above if you wish to see this plot.
axis([350,1000,0,1])
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,10)[2], 'b.',label="10nm" if i==350 else "",ms=1)
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,50)[2], 'r.',label="50nm" if i==350 else "",ms=1)
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,100)[2], 'g.',label="100nm" if i==350 else "",ms=1)
xlabel('Wavelength (nm)')
ylabel('Absorption')
plt.legend(loc="upper left",markerscale=4)
grid()'''