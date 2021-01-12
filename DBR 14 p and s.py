#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 20:28:18 2021

@author: Mark Richmond
"""
from scipy.interpolate import CubicSpline
from scipy.optimize import leastsq
from pylab import *
from numpy.lib.scimath import arcsin
import numpy as np

'''
This program prompts the user to input the variables for a DBR model:
design wavelength, number of periods in the stack, two materials in the stack,
substrate material and medium. It is currently configured to run at two
incident angles (0 and pi/3) for p-polarisation and plots the spectral power
reflectance. Optional line of code for optical resonator is commented below.
'''

def refractive_real_n(material,lambda1):
    text=material.split('\n')
    text[0]+='.txt'
    text ='\n'.join(text)
    w, n, k = loadtxt(text, skiprows=1, unpack=True)
    N = CubicSpline(w, n)
    return (N(lambda1))

def refractive_imaginary_k(material,lambda1):
    text=material.split('\n')
    text[0]+='.txt'
    text ='\n'.join(text)
    w, n, k = loadtxt(text, skiprows=1, unpack=True)
    K = CubicSpline(w, k)
    return (K(lambda1))

incident_angle=0

#Input design parameters of DBR
print('Design wavelength in vacuum(nm) e.g. 633:')
wavelength_des=int(input())
print('Number of periods e.g. 8:')
N = 2*(int(input()))

#Initiaite arrays
Refractive_real=zeros((N+2))  
Refractive_imaginary=zeros((N+2))
R=zeros((N+2),complex)
Layer_thickness=zeros((N+2))
Materials=['' for x in range(N+2)]

print('Medium e.g. Air:')
Materials[0]=input()
Refractive_real[0]=refractive_real_n(Materials[0],wavelength_des)

for i in range(2):
    print('Layer material e.g. Ta2O5 for layer 1 and MgF2 for layer 2:', i+1)
    Materials[i+1]=input()
    Refractive_real[i+1]=refractive_real_n(Materials[i+1],wavelength_des)
    j=1
    while 1+2*j<N:
        Materials[i+1+2*j]=Materials[i+1]
        Refractive_real[i+1+2*j]=Refractive_real[i+1]
        j=j+1
print('Substrate e.g. BK7:')
Materials[N+1]=input()
Refractive_real[N+1]=refractive_real_n(Materials[N+1],wavelength_des)
Layer_thickness=wavelength_des/(4*Refractive_real)

'''
Option to change a central low index layer thickness for optical resonator.
Layer_thickness[8]=230

Function z deals with s- and p-polarisation: select s or p as the 3rd
argument in the call. It also takes incident_angle and wavelength as arguments
and returns the reflection coefficients for a multi-layered stack.
'''

def z(incident_angle,wavelength,polarisation):
    for i in range(N+2):
        Refractive_real[i]=refractive_real_n(Materials[i],wavelength)
        Refractive_imaginary[i]=refractive_imaginary_k(Materials[i],wavelength)
    R=Refractive_real+Refractive_imaginary*(0+1j)
    theta = arcsin((R[0]*sin(incident_angle))/R)
    K = (2*pi/wavelength)*R*cos(theta)
    z=zeros((2,2))
    for i in range(2):
        z[i,i]=1
    for i in range(N):
        T=zeros((2,2),complex)
        if polarisation =='s':
            T[0,0]= 0.5*(sqrt(K[i+1]/K[i])+sqrt(K[i]/K[i+1]))
            T[0,1]= 0.5*(sqrt(K[i+1]/K[i])-sqrt(K[i]/K[i+1]))
        else:
            T[0,0]= 0.5*((R[i+1]/R[i])*sqrt(K[i]/K[i+1])+(R[i]/R[i+1])*sqrt(K[i+1]/K[i]))
            T[0,1]= 0.5*((R[i+1]/R[i])*sqrt(K[i]/K[i+1])-(R[i]/R[i+1])*sqrt(K[i+1]/K[i]))
        T[1,0]=T[0,1]
        T[1,1]=T[0,0]
        P=zeros((2,2),complex)
        P[0,0]=exp((0+1j)*(K[i+1]+0j)*Layer_thickness[i+1])
        P[1,1]=exp((0-1j)*(K[i+1]+0j)*Layer_thickness[i+1])
        x=matmul(P,T)
        z=matmul(x,z)
    Tn=zeros((2,2),complex)
    if polarisation =='s':    
        Tn[0,0]= 0.5*(sqrt(K[N+1]/K[N])+sqrt(K[N]/K[N+1]))
        Tn[0,1]= 0.5*(sqrt(K[N+1]/K[N])-sqrt(K[N]/K[N+1]))
    else:
        Tn[0,0]= 0.5*((R[N+1]/R[N])*sqrt(K[N]/K[N+1])+(R[N]/R[N+1])*sqrt(K[N+1]/K[N]))
        Tn[0,1]= 0.5*((R[N+1]/R[N])*sqrt(K[N]/K[N+1])-(R[N]/R[N+1])*sqrt(K[N+1]/K[N]))        
    Tn[1,0]=Tn[0,1]
    Tn[1,1]=Tn[0,0]
    z=matmul(Tn,z)
    r=-z[1,0]/z[1,1]
    t=z[0,0]+(z[0,1]*r)
    return [r,abs(r)**2]

'''
Plots the reflectance for wavelengths between 350 and 1000nm for
p-polarisation. To change to s-polarisation please change the third argument
in the call to the z function for the relevant plot.
To plot a more detailed graph reduce the number taken as the third argument in
the arange commands below.
'''

thetai=0
for i in arange(350,1000,5):
    plt.plot(i,z(thetai,i,'p')[1], 'b.',label="0" if i==350 else "",ms=2)

thetai=pi/3
for i in arange(350,1000,5):
    plt.plot(i,z(thetai,i,'p')[1], 'r.',label="$\pi$/3" if i==350 else "",ms=2)

xlabel('Wavelength (nm)')
ylabel('Reflectance')
plt.legend(loc="upper right",markerscale=4)
grid()