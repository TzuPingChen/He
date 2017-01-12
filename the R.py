# -*- coding: cp950 -*-
import matplotlib.pyplot as plt  
import mpl_toolkits.mplot3d
import random
import numpy as np
import math
from scipy.spatial.distance import cdist
from visual import*
from visual import*
import sys
from types import*
from time import clock , time
e=math.e
pi=math.pi
k=1.3806448E-23
h=6.62606957E-34

V=1E-27
V4=V**4
Nhe=4
k=1.3806448E-23
m=4*1.66E-27
e=10.22
z=0.2556E-9
L=V**(1.0/3)

def r(a,b):
    r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
    r=r**0.5
    return r
def V(a,b):
    r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
    r=r**0.5
    com=4*e*((z/r)**12-(z/r)**6)
    return com
def R(Nhe):
    a=0
    b=0
    for i in range(Nhe-1):
       c=1-math.cos(2*pi*(i+1)/Nhe)
       a+=c**-6
       b+=c**-3
    com=z**6*a/8/b
    return com**(1.0/6)

for i in range(100):
    print 'Nhe=',i+2,',','R=',R(i+2)
    
