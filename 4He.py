import matplotlib.pyplot as plt
import random

k=1.3806448E-23
m=4*1.66E-27
e=10.22
z=0.2556E-9
L=1E-9
dt=0.0001E-9
N=3

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

x0=[[0,0,0],[0,0,0]]
x1=[[0,0,0],[0,0,0]]
x2=[[0,0,0],[0,0,0]]
x3=[[0,0,0],[0,0,0]]
X=[x0,x1,x2,x3]

w1=[]

for x in X:
    for i in range(3):
        x[0][i]=random.uniform(-L*0.5,L*0.5)
        x[1][i]=random.uniform(-5E-9,5E-9)

for n in range(N):
    w1.append(x[0][1])
    for x in X:
        for i in range(3):
            x[0][i]+=x[1][i]*dt
        while -L<x[0][0]*2<L or -L<x[0][1]*2<L or -L<x[0][2]*2<L:
            for i in range(3):
                if x[0][i]*2>L:
                    x[0][i]=L-x[0][i]
                    x[1][i]=-x[1][i]
                elif x[0][i]*2<-L:
                    x[0][i]=-L-x[0][i]
                    x[1][i]=-x[1][i]
        for y in X:
            if y!=x:
                for j in range(3):
                    x[1][j]-=dV(x,y)*(uni(x,y)[j])*dt
        
plt.plot(w1)
plt.show()

        



