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
import cPickle as pk
import copy

V=1E-27
Nhe=10
Ti=200
MCtr=100000
MCtr2=1
MDtr=10000
dt=1E-25
pi=math.pi
k=1.3806448E-23
beta=1.0/Ti/k
m=6.64E-27
mo_avg=(8*k*Ti*m/pi)**0.5
mo_rms=(3*k*m*Ti)**0.5

def mo_bzm(mo_rms):
    def f(mo):
        c=(2*beta**3/pi/m**3)**0.5*mo**2*math.e**(-beta*mo**2/2/m)/(2*k*Ti*m)**0.5
        return c
    fh=f((2*k*Ti*m)**0.5)
    r=random.random()
    mo=random.random()*mo_rms*3
    while r>(f(mo)/fh):
        mo=random.random()*mo_rms*3
        r=random.random()
    plt.plot(mo,r,'o')
    return mo
for i in range(10):
    mo_bzm(mo_rms)

def f(mo):
    c=(2*beta**3/pi/m**3)**0.5*mo**2*math.e**(-beta*mo**2/2/m)/(2*k*Ti*m)**0.5
    return c
fh=f((2*k*Ti*m)**0.5)
print mo_bzm(mo_rms)
x=np.linspace(0,mo_rms*4,500)
plt.plot(x,f(x)/fh)
plt.show()
