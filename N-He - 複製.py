# -*- coding: cp950 -*-
import matplotlib.pyplot as plt
import random
import math

V=1E-27
VV=[]
PP=[]
    
def run(V,Nhe=4):
    pi=math.pi
    k=1.3806448E-23
    m=4*1.66E-27
    e=10.22
    z=0.2556E-9
    L=V**(1.0/3)
    dt=1E-25   #不可小於mo
    N=10000   #<100000
    cutP=N-1  #<N  尾端計壓力P之範圍
    Ti=20000
    mo=(3*k*Ti*m)**0.5
    p=0
    def R(Nhe):
        a=0
        b=0
        for i in range(Nhe-1):
           c=1-math.cos(2*pi*(i+1)/Nhe)
           a+=c**-6
           b+=c**-3
        com=z**6*a/8/b
        return com**(1.0/6)
    def r(a,b):
        r=(a[0][0]-b[0][0])**2+(a[0][1]-b[0][1])**2+(a[0][2]-b[0][2])**2
        r=r**0.5
        return r
    def V(a,b):
        r=(a[0][0]-b[0][0])**2+(a[0][1]-b[0][1])**2+(a[0][2]-b[0][2])**2
        r=r**0.5
        com=4*e*((z/r)**12-(z/r)**6)
        return com
    def dV(a,b):
        r=(a[0][0]-b[0][0])**2+(a[0][1]-b[0][1])**2+(a[0][2]-b[0][2])**2
        r=r**0.5
        com=24*e*((z**6)*(r**-7)-2*(z**12)*(r**-13))
        return com
    def uni(a,b):
        r=(a[0][0]-b[0][0])**2+(a[0][1]-b[0][1])**2+(a[0][2]-b[0][2])**2
        r=r**0.5
        com=[b[0][0]-a[0][0],b[0][1]-a[0][1],b[0][2]-a[0][2]]
        for i in range(3):
            com[i]=com[i]*1.0/r
        return com
    def T(a):
        com=(a[1][0])**2+(a[1][1])**2+(a[1][2])**2
        com=com/m/2.0
        return com
    
    
    def ml():
        R=random.random()
        th=random.uniform(0,2*math.pi)
        x=R*math.cos(th)
        y=R*math.sin(th)
        def pm():
            if random.random()-0.5>0:
                return 1
            else:
                return -1
        z=pm()*(1-x**2-y**2)**0.5
        return [mo*x,mo*y,mo*z]
    X=[]
    for i in range(Nhe):    #初始動量
        X.append([[0,0,0],ml()])  
    Rhe=R(Nhe)
    for i in range(Nhe):    #初始位置
        X[i][0][0]=Rhe*math.cos(2*pi*i/Nhe)
        X[i][0][1]=Rhe*math.sin(2*pi*i/Nhe)                                                                                                            
    
    AllN=[]
    for i in range(Nhe):
        AllN.append([[],[],[]])  #每粒資訊:[[x...],[y...],[z...]]
    pp=[]
    AllU=[]
    AllT=[]
    AllV=[]
    p=0
    def rec():   #每步紀錄
        for i in range(Nhe):
            for j in range(3):
                AllN[i][j].append(X[i][0][j])
        AlV=0
        for i in range(Nhe):
            for j in range(Nhe-1-i):
                AlV+=V(X[i],X[Nhe-1-j])
        
        AlT=0
        for i in range(Nhe):
            AlT+=T(X[i])
        AllV.append(AlV)
        AllT.append(AlT)
        pp.append(p)
    
    rec()
    for n in range(N):
        for x in X:
            for y in X:
                if y!=x:
                    for j in range(3):
                        x[1][j]+=dV(x,y)*(uni(x,y)[j])*dt
        for x in X:
            for i in range(3):
                x[0][i]+=x[1][i]*dt/m
        for x in X:
            for i in range(3):
                if x[0][i]>L*0.5:
                    x[0][i]=L-x[0][i]
                    x[1][i]=-x[1][i]
                    p+=(x[1][i]**2)**0.5
                elif x[0][i]<-L*0.5:
                    x[0][i]=-x[0][i]-L
                    x[1][i]=-x[1][i]
                    p+=(x[1][i]**2)**0.5
                else:
                    p+=0
        rec()
    for i in range(N-1):
        AllU.append(AllT[i]+AllV[i+1])

    plt.subplot(331)
    for i in range(Nhe):
        plt.plot(AllN[i][0],AllN[i][1])
        
    plt.subplot(332)
    for i in range(Nhe):
        plt.plot(AllN[i][0],AllN[i][2])

    plt.subplot(333)
    plt.plot(pp)

    plt.subplot(334)
    plt.plot(AllU)

    plt.subplot(338)
    plt.plot(AllT)

    plt.subplot(339)
    plt.plot(AllV)
    
    print "V=",L**3,"(L=",L,")"
    print "Ti=",Ti,"(per T=",mo**2/2/m,")"
    print "T+V=U=",AllU[0]
    print 'P=',2*(pp[N]-pp[N-cutP])/L**2/6/cutP/dt
    plt.show()

run(V)

