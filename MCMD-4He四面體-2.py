# -*- coding: cp950 -*-
import matplotlib.pyplot as plt  
import mpl_toolkits.mplot3d
import random
import math
import copy
import numpy as np

def go(V,Ti):
    Nhe=4    #粒子數
    k=1.3806448E-23
    MCtr=10000   #MC1步數
    MCtr2=1       #MC2步數
    MDtr=10000    #MD 步數
    dt=1E-25
    pi=math.pi
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
        return mo,r
    def f(mo):
        c=(2*beta**3/pi/m**3)**0.5*mo**2*math.e**(-beta*mo**2/2/m)/(2*k*Ti*m)**0.5
        return c
    fh=f((2*k*Ti*m)**0.5)
    mo=mo_bzm(mo_rms) #x初始動量選擇(mo_avg,mo_rms,mo_bzm(mo_rms))
    e=10.22
    z=0.2556E-9
    beta=1.0/k/Ti

    L=V**(1.0/3)
    Lo=2**(-1.0/3)*z
    X_=[]
    #粒子初始位置設定(判斷是否壓縮四面體)
    if L>Lo:
        X_.append([-Lo*0.5,-Lo*0.5,-Lo*0.5])
        X_.append([Lo*0.5,Lo*0.5,-Lo*0.5])
        X_.append([-Lo*0.5,Lo*0.5,Lo*0.5])
        X_.append([Lo*0.5,-Lo*0.5,Lo*0.5])
    else:
        X_.append([-L*0.5,-L*0.5,-L*0.5])
        X_.append([L*0.5,L*0.5,-L*0.5])
        X_.append([-L*0.5,L*0.5,L*0.5])
        X_.append([L*0.5,-L*0.5,L*0.5])
    oldX=copy.copy(X_)
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
    print 'Vo=',Lo**3,'Lo=',Lo
    print 'MCtry=',MCtr
    ax = plt.subplot(111,projection='3d',title='Monte Carlo')
    ax.plot(oxx,oyy,ozz,'o',c='red',label='start place')
    ax.plot(xx,yy,zz,'o',c='blue',label='equilibrium place\n(for MD\'s start place)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.legend()
    plt.show()

    def run(X_,V,Nhe,L,Ti,pi,k,m,e,z,dt,MDtr,mo):
        cutP=MDtr-1  #<N  尾端計壓力P之範圍
        p=0
        def r(a,b):
            r=(a[0][0]-b[0][0])**2+(a[0][1]-b[0][1])**2+(a[0][2]-b[0][2])**2
            r=r**0.5
            return r
        def Vf(a,b):
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
        for i in range(Nhe):    #初始動量
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
                    AlV+=Vf(X[i],X[Nhe-1-j])
            
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
        print "Uk=",AllT[0]
        print "Vtot=",AllV[0],',(Vtot_min=',-6.0*e,')'
        print "Utot=",AllU[0]
        print 'P=',2*(pp[MDtr]-pp[MDtr-cutP])/L**2/6/cutP/dt
        
    run(X_,V,Nhe,L,Ti,pi,k,m,e,z,dt,MDtr,mo)
    x=np.linspace(0,mo_rms*4,500)
    plt.subplot(337)
    plt.plot(x,f(x)/fh)
    plt.show()
k=1.3806448E-23
z=0.2556E-9
V=1.01*(2**(-1.0/3)*z)**3
Ti=1.0/k/1.0
go(V,Ti)
