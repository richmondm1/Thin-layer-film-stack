#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:31:42 2021

@author: Mark Richmond
"""
from scipy.interpolate import CubicSpline
from pylab import *
import cmath

'''
This program plots the reflectance spectrum for a single layer film on a
substrate in a given medium at varied angles of incidence. It is set up for a
MgF2 layer on a substrate of BK7 in air.
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
reflectance as calculated from the appropriate matrix elements for a varying
angle of incidence. Equations for s-polarization are used. Single layer
between medium and substrate.
'''

def z(wavelength,layer_thickness):
    Refractive_real = array([N0(wavelength),N1(wavelength),N2(wavelength)])
    Refractive_imaginary = array([K0(wavelength),K1(wavelength),K2(wavelength)])*(0+1j)
    R = Refractive_real+Refractive_imaginary
    theta=arcsin((R[0]*sin(incident_angle))/R)
#K is now the z-component of the wavevector in each layer.
    K=((2*pi)/wavelength)*R*cos(theta)
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
Plots the reflectance for wavelengths between 350 and 1000nm for three
different angles: 0, pi/6, pi/3. To plot a more detailed graph reduce the
number taken as the third argument in the arange commands below.
'''

axis([350,1000,0,0.18])

incident_angle=0
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,90), 'b.',label="0" if i==350 else "",ms=1)
    
incident_angle=pi/6
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,90), 'r.',label="$\pi$/6" if i==350 else "",ms=1)
    
incident_angle=pi/3
for i in arange(350,1000,0.5):
    plt.plot(i,z(i,90), 'g.',label="$\pi$/3" if i==350 else "",ms=1)

xlabel('Wavelength (nm)')
ylabel('Reflectance')
plt.legend(loc="upper left",markerscale=4)
grid()
