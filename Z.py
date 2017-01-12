import random
import numpy as np
import math
import matplotlib.pyplot as plt
V=1E-27
Nhe=100
Ti=200
pi=math.pi
k=1.3806448E-23
h=6.626E-34
e=10.22
z=0.2556E-9
m=6.64E-27
beta=1.0/k/Ti
Ztr=10000

def Zqa(Ti,V,Nhe,k,e,Ztr):
    def je(x):
        c=1
        for i in range(x):
           c=c*(i+1)
        return c
    X=[]
    L=V**(1.0/3)
    beta=1.0/k/Ti
    C_=[]
    Cqa=0
    exp=math.e
    amax=Nhe*(Nhe-1)*0.5*e*beta
    for abc in range(Ztr):
        X=[]
        for i in range(Nhe):  
            X.append([random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5)])
        def Vf(a,b):
            r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
            r=r**0.5
            com=4*e*((z/r)**12-(z/r)**6)
            return com
        Cvq=0
        for i in range(Nhe):
            for j in range(Nhe-1-i):
                Cvq+=Vf(X[i],X[Nhe-1-j])
        C_.append(-1.0*beta*Cvq)
    Cvqmax=max(C_)
    for i in range(Ztr):
        Cqa+=exp**(C_[i]-Cvqmax)
    return Nhe*math.log(V)+Cvqmax+math.log(Cqa/Ztr)


print Zqa(Ti,V,Nhe,k,e,Ztr)
print -beta*(p**2.0/m+V)-3.0*N/2*math.log(2*pi*m*k*Ti)

"""
def P(p,V,beta,m,pi,k,Ti,Zq):
    up=-beta*(p**2.0/m+V)-3.0*N/2*math.log(2*pi*m*k*Ti)-lZq
    return exp**(up)
"""
