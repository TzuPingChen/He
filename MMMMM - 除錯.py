# -*- coding: cp950 -*-
import matplotlib.pyplot as plt  
import mpl_toolkits.mplot3d
import random
import math
import copy
def U(Ti,V,Nhe,MCtr):  
    pi=math.pi
    exp=math.e
    k=1.3806448E-23
    e=10.22
    z=0.2556E-9
    L=V**(1.0/3)
    XX=[] #共MCtr個X系統
    XE=[] #每系統對應能量
    Ee=[]
    def Vtot(X):  #系統總位能
        def Vf(a,b):
            r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
            r=r**0.5
            com=4*e*((z/r)**12-(z/r)**6)
            return com
        com=0
        for i in range(Nhe):
                for j in range(Nhe-1-i):
                    com+=Vf(X[i],X[Nhe-1-j])
        return com
    for i in range(MCtr):
        X=[]
        for j in range(Nhe):
            X.append([random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5)])
        XE.append(Vtot(X))
        XX.append(X)
    Emin=min(XE)
    for i in range(MCtr):
        Ee.append(exp**(-beta*(XE[i]-Emin)))
    up=0
    for i in range(MCtr):
        up+=XE[i]*Ee[i]
    down=sum(Ee)
    Eavg=up*1.0/down/Nhe
    print 'Emin=',Emin/Nhe
    print 'Eavg=',Eavg
    #print 'XX=',XX
    #print 'XE=',XE
    #print 'Ee=',Ee
    #print 'up(up+=XE[i]*Ee[i])=',up
    #print 'down(sumEe)=',down

V=1E-27
Nhe=8
beta=0.1
k=1.3806448E-23
Ti=1.0/beta/k
MCtr=10000
U(Ti,V,Nhe,MCtr)

'''
for i in range(10):
    Ti=1.0/beta/k
    print 'beta=',U(Ti,V,Nhe,MCtr)[0]
    print 'Emin=',U(Ti,V,Nhe,MCtr)[1]
    print 'Eavg=',U(Ti,V,Nhe,MCtr)[2]
    plt.subplot(111)
    plt.plot(U(Ti,V,Nhe,MCtr)[0],U(Ti,V,Nhe,MCtr)[2],'+')
    plt.subplot(111)
    plt.plot(U(Ti,V,Nhe,MCtr)[0],U(Ti,V,Nhe,MCtr)[1],'.')
    beta+=0.1
plt.show()
    
'''
