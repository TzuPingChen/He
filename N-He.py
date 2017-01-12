import matplotlib.pyplot as plt
import random
import math

k=1.3806448E-23
m=4*1.66E-27
e=10.22
z=0.2556E-9
V=8E-27
L=V**(1.0/3)
s=1E-9
dt=1E-27
N=50000
mo=1E-11
p=0

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

def hl():
    com=random.uniform(s*0.25,s*0.5)
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

x0=[[0,0,0],[0,0,0]]
x1=[[0,0,0],[0,0,0]]
x2=[[0,0,0],[0,0,0]]
x3=[[0,0,0],[0,0,0]]
X=[x0,x1,x2,x3]

x0[0]=[hl(),hl(),hl()]
x1[0]=[-hl(),-hl(),hl()]
x2[0]=[-hl(),hl(),-hl()]
x3[0]=[hl(),-hl(),-hl()]
for x in X:
    x[1]=ml()

ax=[]
ay=[]
az=[]
aT=[]
bx=[]
by=[]
bz=[]
bT=[]
cx=[]
cy=[]
cz=[]
cT=[]
dx=[]
dy=[]
dz=[]
dT=[]
abV=[]
acV=[]
adV=[]
bcV=[]
bdV=[]
cdV=[]
uu=[]
pp=[]

def rec():
    ax.append(x0[0][0])
    ay.append(x0[0][1])
    az.append(x0[0][2])
    bx.append(x1[0][0])
    by.append(x1[0][1])
    bz.append(x1[0][2])
    cx.append(x2[0][0])
    cy.append(x2[0][1])
    cz.append(x2[0][2])
    dx.append(x3[0][0])
    dy.append(x3[0][1])
    dz.append(x3[0][2])
    abV.append(V(x0,x1))
    acV.append(V(x0,x2))
    adV.append(V(x0,x3))
    bcV.append(V(x1,x2))
    bdV.append(V(x1,x3))
    cdV.append(V(x2,x3))
    uu.append(T(x0)+T(x1)+T(x2)+T(x3)+V(x0,x1)+V(x0,x2)+V(x0,x3)+V(x1,x2)+V(x1,x3)+V(x2,x3)) 
    aT.append(T(x0))
    bT.append(T(x1))
    cT.append(T(x2))
    dT.append(T(x3))
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
                x[1][i]=-x[1][i]
                p+=(x[1][i]**2)**0.5
            if x[0][i]<-L*0.5:
                x[1][i]=-x[1][i]
                p+=(x[1][i]**2)**0.5
    rec()

print 'V=',L**3
print 'P=',2*pp[N]/L**2/6/N/dt


plt.subplot(331)
plt.plot(ax,ay)
plt.plot(bx,by)
plt.plot(cx,cy)
plt.plot(dx,dy)

plt.subplot(333)
plt.plot(aT)

plt.subplot(334)
plt.plot(bT)

plt.subplot(335)
plt.plot(cT)

plt.subplot(336)
plt.plot(dT)

plt.subplot(339)
plt.plot(uu)

plt.subplot(338)
plt.plot(pp)

plt.show()
