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
exp=math.e
pi=math.pi
k=1.3806448E-23
h=6.62606957E-34

V=1E-27
V4=V**4
N=4     #����
k=1.3806448E-23
m=4*1.66E-27
e=10.22
z=0.2556E-9
L=V**(1.0/3)
Ee=0
End=1000
dE=1
tr=100000
omq=[[],[]]

def r(a,b):
    r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
    r=r**0.5
    return r

def Vup(a,b):
    r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
    r=r**0.5
    com=4*e*((z/r)**12-(z/r)**6)
    return com+e

def V(a,b):
    r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
    r=r**0.5
    com=4*e*((z/r)**12-(z/r)**6)
    return com
X=[]
for i in range(N):
    X.append([0,0,0])
    
his=zeros(End/dE,int) 
for n in range(tr):
    for x in X:
        for i in range(3):
            x[i]=L*(random.random()-0.5)
    sumV=0
    for i in range(N):
            for j in range(N-1-i):
                sumV+=Vup(X[i],X[N-1-j])        
    if sumV<End:
        his[sumV//dE]+=1    
def en(N):
    com=N*(1-N)*0.5
    return com*e  

for i in range(int(End/dE)):
    omq[0].append((0.5+i)*dE+en(N))
    if his[i]==0:
        his[i]+=1
    else:
        omq[1].append(his[i])
    
    
def omp(E):
    def ng(x):
        a=1
        for i in range(x):
            a=a*(i+1)
        return a*1.0
    com=3*N*pi**(3*N/2.0)*m*(2*m*E)**((3*N-2)/2.0)/ng(3*N/2)
    return com

    
def omega_(E,dE,N):
    En=E-en(N)
    s=0
    for i in range(int(En/dE)):
        s+=omq[1][int(En/dE)-i-1]*omp((0.5+i)*dE)*dE
    return s

"""
GOp=[]
GOq=[]
E=60
for i in range(int(E/dE)):
    GOp.append(omp((0.5+i)*dE))
    GOq.append(omq[1][int(E/dE)-i-1])


plt.subplot(221)
plt.plot(GOp)

plt.subplot(222)
plt.plot(GOq)
plt.show()

"""
print en(N)
print omega_(-50,dE,N)
print omega_(-50.2,dE,N)

"""

EE=omq[0]
OO=[]
dOO=[]
beta=[]
TT=[]
for i in range(int(End/dE)):
    x=(1+i)*dE+en(N)
    OO.append(omega_(x,dE,N))

dOO.append(0)
for i in range(int(End/dE)-1):
    dOO.append((OO[i+1]-OO[i])*1.0/dE)
beta.append(0)
for i in range(int(End/dE)-1):
    if OO[i+1]==0:
        beta.append(0)
    else:
        beta.append(dOO[i]/OO[i+1])
TT.append(0)
for i in range(int(End/dE)-1):
    if beta[i]==0:
        TT.append(0)
    else:
        TT.append(1.0/beta[i]/k)

OT=[]
def Ot(E,N):
    com=2.0*E/3/N/k
    return com
for i in range(len(EE)):
    OT.append(Ot(EE[i],N))

plt.subplot(331)
plt.plot(EE,OO,'red')

plt.subplot(332)
plt.plot(EE,dOO)

plt.subplot(333)
plt.plot(EE,beta)

plt.subplot(334)
plt.plot(EE,TT)

plt.subplot(337)
plt.plot(EE,OT)

plt.subplot(336)
plt.plot(his)

plt.show()

"""

"""
def beta(E):
    com=omega_(E+dE*0.5,End,dE,tr)-omega_(E-dE*0.5,End,dE,tr)
    com=com/dE
    return com

print beta(3.3)
"""

