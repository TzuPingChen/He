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
k=1.3806448E-23
h=6.62606957E-34
def omega_(E,End=1000,dE=0.01,tr=1000):   #omega(E)=(1/h^12/24)*omega_(E)  ,but(1/h^12/24) too large
    pi=math.pi
    
    V=8E-27
    V4=V**4
    N=4
    k=1.3806448E-23
    m=4*1.66E-27
    e=10.22
    z=0.2556E-9
    L=V**(1.0/3)
    Ee=0
    def r(a,b):
        r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
        r=r**0.5
        return r
    def V(a,b):
        r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
        r=r**0.5
        com=4*e*((z/r)**12-(z/r)**6)
        return com+e
        
    x0=[0,0,0]
    x1=[0,0,0]
    x2=[0,0,0]
    x3=[0,0,0]
    X=[x0,x1,x2,x3]
    qv=[]
        
    his=zeros(End/dE,int) ##0~999
        
    for n in range(tr):
        for x in X:
            for i in range(3):
                x[i]=L*(random.random()-0.5)
        sumV=V(x0,x1)+V(x0,x2)+V(x0,x3)+V(x1,x2)+V(x1,x3)+V(x2,x3)
        if sumV<End:
            his[sumV//dE]+=1
    
    
    omq=[[],[]]
    for i in range(int(E/dE)):
        omq[0].append((0.5+i)*dE)
        omq[1].append(his[i]*V4*1.0/tr)
    
    
    def omp(E):
        com=3*N*pi**(3*N/2.0)*m*(2*m*E)**((3*N-2)/2.0)/720.0
        return com
    com=0
    


    for i in range(int(E/dE)):
        com+=omq[1][int(E/dE)-i-1]*omp((0.5+i)*dE)*dE
    return com


EE=[]
OO=[]
End=500
dE=1
dOO=[]
beta=[]
TT=[]

for i in range(int(End/dE)):
    x=(1+i)*dE
    EE.append(x)
    OO.append(omega_(x,End))

dOO.append(0)
for i in range(int(End/dE)-1):
    dOO.append((OO[i+1]-OO[i])*1.0/dE)
beta.append(dOO[0]/OO[0+1])
for i in range(int(End/dE)-1):
    beta.append(dOO[i]/OO[i+1])
TT.append(1000)
for i in range(int(End/dE)-1):
    if beta[i]==0:
        TT.append(1000)
    else:
        TT.append(1.0/beta[i]/k)


plt.subplot(221)
plt.plot(EE,OO)

plt.subplot(222)
plt.plot(EE,dOO)

plt.subplot(223)
plt.plot(EE,beta)

plt.subplot(224)
plt.plot(EE,TT)
plt.show()

"""
def beta(E):
    com=omega_(E+dE*0.5,End,dE,tr)-omega_(E-dE*0.5,End,dE,tr)
    com=com/dE
    return com

print beta(3.3)
"""

