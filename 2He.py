import matplotlib.pyplot as plt
import random

k=1.3806448E-23
m=4*1.66E-27
e=10.22
z=0.2556E-9
L=2E-9
s=0.5E-9
dt=1E-27
N=10000
mo=1E-10
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
x0=[[0,0,0],[0,0,0]]
x1=[[0,0,0],[0,0,0]]
X=[x0,x1]

ax=[]
ay=[]
az=[]
abV=[]
aT=[]
bx=[]
by=[]
bz=[]
bT=[]
rr=[]
ddV=[]
uu=[]

for i in range(3):
    x0[0][i]=random.uniform(s*0.25,s*0.5)
    x1[0][i]=random.uniform(-s*0.25,-s*0.5)
    x0[1][i]=random.uniform(-mo,mo)
    x1[1][i]=random.uniform(-mo,mo)
    

def rec():
    rr.append(r(x0,x1))
    ax.append(x0[0][0])
    ay.append(x0[0][1])
    az.append(x0[0][2])
    bx.append(x1[0][0])
    by.append(x1[0][1])
    bz.append(x0[0][2])
    abV.append(V(x0,x1))
    ddV.append(dV(x0,x1))
    uu.append(T(x0)+T(x1)+V(x0,x1)) 
    aT.append(T(x0))
    bT.append(T(x1))

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
            if x[0][i]<-L*0.5:
                x[1][i]=-x[1][i]
    rec()



plt.subplot(331)
plt.plot(ax,ay)
plt.plot(bx,by)

plt.subplot(333)
plt.plot(rr,label="distance")

plt.subplot(334)
plt.plot(aT,label="aT")

plt.subplot(335)
plt.plot(bT,label="bT")

plt.subplot(337)
plt.plot(abV,label="abV")

plt.subplot(339)
plt.plot(uu,label="U")

plt.show()


        



