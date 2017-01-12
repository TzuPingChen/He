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
e=10.22
def en(N):
    com=N*(1-N)*0.5
    return com*e  

V=1000E-27
V4=V**4
N=4     #偶數
k=1.3806448E-23
m=4*1.66E-27
z=0.2556E-9
L=V**(1.0/3)
Ee=0
End=0.5E3-en(N)
dE=0.1
tr=100000
omq=[[],[]]
End-=en(N)
Ti=200
beta=1.0/Ti/k


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
    if sumV<End-1:
        his[sumV//dE]+=1    

for i in range(int(End/dE)):
    omq[0].append((0.5+i)*dE+en(N))
    if his[i]==0:
        his[i]+=0
        omq[1].append(his[i]*1E-6)
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
EE=[]
OO=[]
PP=[]
PPP=[]
EEE=[]
ppp=[]
TRY=int(End/dE)
for i in range(TRY):
    EE.append((en(N)+1)/dE+i*dE)
    OO.append(omega_((en(N)+1)/dE+i*dE,dE,N))
    if omega_((en(N)+1)/dE+i*dE,dE,N)==0:
        PP.append(0)
    else:
        PP.append(math.log(omega_((en(N)+1)/dE+i*dE,dE,N))+(-beta*((en(N)+1)/dE+i*dE)))
for j in range(int(End/dE)):
    if PP[j]!=0:
        EEE.append(EE[j])
        ppp.append(PP[j])
PPmax=max(ppp)
for j in range(len(ppp)):
    PPP.append(exp**(ppp[j]-PPmax))
 
plt.subplot(221)
plt.plot(EE,OO)
plt.subplot(222)
plt.plot(omq[1])
plt.subplot(223)
plt.plot(EE,PP)
plt.subplot(224)
plt.plot(EEE,PPP)
plt.show()


























































