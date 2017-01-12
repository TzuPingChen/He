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

def omega(E):
    
    pi=math.pi
    h=6.62606957E-34
    
    V=8E-27
    V4=V**4
    N=4
    k=1.3806448E-23
    m=4*1.66E-27
    e=10.22
    z=0.2556E-9
    L=V**(1.0/3)
    tr=50000
    dE=1
    Ee=0
    End=1000  #需整數
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
    
    his=zeros(End,int)
    
    for n in range(tr):
        for x in X:
            for i in range(3):
                x[i]=L*(random.random()-0.5)
        sumV=V(x0,x1)+V(x0,x2)+V(x0,x3)+V(x1,x2)+V(x1,x3)+V(x2,x3)
        if sumV<End:
            his[sumV//dE]+=1
    
    omq=[[],[]]
    for i in range(End):
        omq[0].append((0.5+i)*dE)
        omq[1].append(his[i]*V4/tr)
    
    def omp(E):
        com=3*N*pi**(3*N/2.0)*m*(2*m*E)**((3*N-2)/2.0)/720.0
        return com
    com=0
    def li(E):
        if E%dE==0:
            return E//dE
        else:
            return E//dE+1
    for i in range(len(his)):
        com+=omq[1][len(his)-i]*omp((0.5+i)*dE)*dE
    return com/h**12/24
        
print omega(1.5)
