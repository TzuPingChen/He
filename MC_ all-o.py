# -*- coding: cp950 -*-
import matplotlib.pyplot as plt  
import mpl_toolkits.mplot3d
import random
import math
import copy

V=1E-27
Nhe=10
pi=math.pi
k=1.3806448E-23
e=10.22
z=0.2556E-9
L=V**(1.0/3)
MCtr=1000
MCtr2=1
Ti=200
beta=1.0/k/Ti
X=[]

for i in range(Nhe):  
    X.append([random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5)])
oldX=copy.copy(X)
def Vtot(n):  #系統總位能
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
AllV=[]
for i in range(MCtr):
    n=random.randint(0,Nhe-1)
    for j in range(MCtr2):
        old=copy.copy(X[n])
        oldV=copy.copy(Vtot(n))
        X[n]=[random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5),random.uniform(-L*0.5,L*0.5)]
        if Vtot(n)>oldV:
            r=random.random()
            if r>math.e**(-beta*(Vtot(n)-oldV)):
                X[n]=old
    def Vf(a,b):
        r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
        r=r**0.5
        com=4*e*((z/r)**12-(z/r)**6)
        return com
    AlV=0
    for i in range(Nhe):
        for j in range(Nhe-1-i):
            AlV+=Vf(X[i],X[Nhe-1-j])
    AllV.append(AlV)

oxx,oyy,ozz,xx,yy,zz=[],[],[],[],[],[]
for i in range(Nhe):
    xx.append(X[i][0])
    yy.append(X[i][1])
    zz.append(X[i][2])
    oxx.append(oldX[i][0])
    oyy.append(oldX[i][1])
    ozz.append(oldX[i][2])
print 'T=',Ti
print 'V=',V,'(L=',L,')'
print 'N=',Nhe
print 'D=N/V=',Nhe*1.0/V
print 'tryN=',MCtr
ax = plt.subplot(221,projection='3d')
ax.plot(oxx,oyy,ozz,'o',c='red')
ax.plot(xx,yy,zz,'o',c='blue')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.subplot(222)
plt.plot(AllV)
plt.show()
