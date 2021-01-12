#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 11:28:58 2021

@author: Mark Richmond
"""
from scipy.interpolate import CubicSpline
from scipy.optimize import leastsq
from pylab import *
from numpy.lib.scimath import arcsin

'''
This program prompts the user to input the variables for a DBR model:
design wavelength, number of periods in the stack, two materials in the stack,
substrate material and medium. It is currently configured to run at zero
incident angle and prints amplitude and power reflection coefficients. It was
used to look at the variation in reflectance with number of periods and also
to study the phase of the reflected wave.
'''

def refractive_real(material):
    text=material.split('\n')
    text[0]+='.txt'
    text ='\n'.join(text)
    w, n, k = loadtxt(text, skiprows=1, unpack=True)
    N = CubicSpline(w, n)
    return (N(wavelength))

def refractive_imaginary(material):
    text=material.split('\n')
    text[0]+='.txt'
    text ='\n'.join(text)
    w, n, k = loadtxt(text, skiprows=1, unpack=True)
    K = CubicSpline(w, k)
    return (K(wavelength))

incident_angle=0
print('Design wavelength (nm) e.g. 633:')
wavelength=int(input())
print('Number of periods e.g. 8:')
N = 2*(int(input()))
Refractive_real=zeros((N+2))  
Refractive_imaginary=zeros((N+2))
R=zeros((N+2),complex)
Layer_thickness = zeros((N+2))
print('Medium e.g. Air:')
g=input()
Refractive_real[0]=refractive_real(g)
Refractive_imaginary[0]=refractive_imaginary(g)
Layer_thickness[0]=wavelength/(4*Refractive_real[0])
for i in range(2):
    print('Layer material e.g. Ta2O5 for layer 1 and MgF2 for layer 2:', i+1)
    y=input()
    j=0
    while 1+2*j<N:
        Refractive_real[i+1+2*j]=refractive_real(y)
        Refractive_imaginary[i+1+2*j]=refractive_imaginary(y)
        Layer_thickness[i+1+2*j]=wavelength/(4*Refractive_real[i+1+2*j])
        j=j+1
print('Substrate e.g. BK7:')
h=input()
Refractive_real[N+1]=refractive_real(h)
Refractive_imaginary[N+1]=refractive_imaginary(h)
R=Refractive_real+(Refractive_imaginary*(0+1j))

#for electric polarised reflection and transmission
def z(incident_angle):
    K_i = (2*pi/wavelength)*R[0]*cos(incident_angle)
    theta_j = arcsin((R[0]*sin(incident_angle))/R[1])
    K_j = (2*pi/wavelength)*R[1]*cos(theta_j)
    z=zeros((2,2))
    for i in range(2):
        z[i,i]=1
    for i in range(N):
        T=zeros((2,2),complex)
        T[0,0]= 0.5*(sqrt(K_j/K_i)+sqrt(K_i/K_j))
        T[0,1]= 0.5*(sqrt(K_j/K_i)-sqrt(K_i/K_j))
        T[1,0]=T[0,1]
        T[1,1]=T[0,0]
        P=zeros((2,2),complex)
        P[0,0]=exp((0+1j)*(K_j+0j)*Layer_thickness[i+1])
        P[1,1]=exp((0-1j)*(K_j+0j)*Layer_thickness[i+1])
        x=matmul(P,T)
        z=matmul(x,z)
        K_i=K_j
        incident_angle = arcsin((R[0]*sin(incident_angle))/R[i+2])
        K_j = (2*pi/wavelength)*R[i+2]*cos(theta_j)
    Tn=zeros((2,2),complex)
    Tn[0,0]= 0.5*(sqrt(K_j/K_i)+sqrt(K_i/K_j))
    Tn[0,1]= 0.5*(sqrt(K_j/K_i)-sqrt(K_i/K_j))
    Tn[1,0]=Tn[0,1]
    Tn[1,1]=Tn[0,0]
    z=matmul(Tn,z)
    r=-z[1,0]/z[1,1]
    t=z[0,0]+(z[0,1]*r)
    print ('printing [r, R]')
    return r,abs(r)**2

print (z(incident_angle))