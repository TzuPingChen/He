# -*- coding: cp950 -*-
import matplotlib.pyplot as plt  
import mpl_toolkits.mplot3d
import random
import math
import copy
import numpy as np

V=1E-20
Nhe=10
Ti=200
MCtr=100000
MCtr2=1
MDtr=10000
Ztr=10000
dt=1E-25
pi=math.pi
k=1.3806448E-23
beta=1.0/Ti/k
m=6.64E-27
exp=math.e
e=10.22
z=0.2556E-9
L=V**(1.0/3)
beta=1.0/k/Ti
X_=[]
for i in range(Nhe):  
    X_.append([random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5)])
oldX=copy.copy(X_)
def Vtot(n):  #第n粒所感位能
    def V(a,b):
        r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
        r=r**0.5
        com=4*e*((z/r)**12-(z/r)**6)
        return com
    com=0
    for i in range(Nhe):
        if n!=i:
            com+=V(X_[n],X_[i])
    return com

for i in range(MCtr):
    n=random.randint(0,Nhe-1)
    for j in range(MCtr2):
        old=copy.copy(X_[n])
        oldV=copy.copy(Vtot(n))
        X_[n]=[random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5)]
        if Vtot(n)>oldV:
            r=random.random()
            if r>math.e**(-beta*(Vtot(n)-oldV)):
                X_[n]=old
oxx,oyy,ozz,xx,yy,zz=[],[],[],[],[],[]
for i in range(Nhe):
    xx.append(X_[i][0])
    yy.append(X_[i][1])
    zz.append(X_[i][2])
    oxx.append(oldX[i][0])
    oyy.append(oldX[i][1])
    ozz.append(oldX[i][2])
print 'T=',Ti
print 'V=',V,'(L=',L,')'
print 'N=',Nhe
print 'D=N/V=',Nhe*1.0/V
print 'MCtry=',MCtr
"""
ax = plt.subplot(111,projection='3d')
ax.plot(oxx,oyy,ozz,'o',c='red')
ax.plot(xx,yy,zz,'o',c='blue')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
"""
def lnZqa(Ti,V,Nhe,k,e,Ztr,exp): #Zqa(T,V,N)  run外計算
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
lnZqa_=lnZqa(Ti,V,Nhe,k,e,Ztr,exp)
print lnZqa_
def run(X_,V,Nhe,L,Ti,pi,k,m,e,z,dt,MDtr,mo,lnZqa_,exp):
    cutP=MDtr-1  #<N  尾端計壓力P之範圍
    p=0    
    def r(a,b):
        r=(a[0][0]-b[0][0])**2+(a[0][1]-b[0][1])**2+(a[0][2]-b[0][2])**2
        r=r**0.5
        return r
    def Vf(a,b):  #待改其他V→Vf
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

    def P(mo,n,beta,m,Nhe,V,k,e,Ztr,Ti,lnZqa_,exp): #第n粒含動能mo之機率
        def Vtot(n):  #第n粒所感位能
            def Vf(a,b):
                r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
                r=r**0.5
                com=4*e*((z/r)**12-(z/r)**6)
                return com
            com=0
            for i in range(Nhe):
                if n!=i:
                    com+=Vf(X_[n],X_[i])
            return com
        zBHup=-beta*(mo**2/m+Vtot(n))
        zpup=-3.0*Nhe/2*math.log(2*pi*m*k*Ti)        
        zqup=lnZqa_
        return (zBHup+zpup-zqup)
    return P(mo,n,beta,m,Nhe,V,k,e,Ztr,Ti,lnZqa_,exp)
mo=0
mozeP=run(X_,V,Nhe,L,Ti,pi,k,m,e,z,dt,MDtr,0,lnZqa_,exp)
for i in range(150):
    print exp**(run(X_,V,Nhe,L,Ti,pi,k,m,e,z,dt,MDtr,mo,lnZqa_,exp)-mozeP)
    mo+=0.001*(i+1)*(8.0*k*Ti*m/pi)**0.5


"""

    def ml():  #初始動量
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
        m0=mo_bzm(mo_rms)
        plt.subplot(337)
        plt.plot(m0[0],m0[1],'o')
        return [m0[0]*x,m0[0]*y,m0[0]*z]
    X=[]
    for i in range(Nhe):
        X.append([X_[i],ml()])                                                                                                            
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
    for n in range(MDtr):
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
    for i in range(MDtr-1):
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

    plt.subplot(335)
    plt.plot(AllT)

    plt.subplot(336)
    plt.plot(AllV)
    print
    print "V=",L**3,"(L=",L,")"
    print "T=",Ti
    print "T+V=U=",AllU[0]
    print 'P=',2*(pp[MDtr]-pp[MDtr-cutP])/L**2/6/cutP/dt
    
run(X_,V,Nhe,L,Ti,pi,k,m,e,z,dt,MDtr,mo)
x=np.linspace(0,mo_rms*4,500)
plt.subplot(337)
plt.plot(x,f(x)/fh)
plt.show()
"""
