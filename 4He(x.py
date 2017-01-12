import matplotlib.pyplot as plt
import random

k=1.3806448E-23
m=4*1.66E-27
e=10.22
z=0.2556E-9
N=3


def V(a,b):
    r=(a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2
    r=r**0.5
    if r>0.22:
        com=4*e*((z/r)**12-(z/r)**6)
    else:
        return 0
x0=[]
x1=[]
x2=[]
x3=[]
X=[x0,x1,x2,x3]
for x in X:
    for i in range(N):
        for j in range(N):
            for k in range(N):
                x.append([i*1.0E-9/N,j*1.0E-9/N,k*1.0E-9/N])
for a in range(N**3):
    for b in range(N**3):
        for c in range(N**3):
            for d in range(N**3):
                Ev=0
                Ev+=V(x0[a],x1[b])
                Ev+=V(x0[a],x2[c])
                Ev+=V(x0[a],x3[d])
                Ev+=V(x1[b],x2[c])
                Ev+=V(x1[b],x3[d])
                Ev+=V(x2[c],x3[d])
                plt.plot(random.random(),Ev)
plt.show()





