# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 14:27:25 2014

@author: Sollphin
"""
import numpy as np
import random as rd
import matplotlib.pyplot as plt
import scipy.constants as scc
import math
import mpl_toolkits.mplot3d

a = 2
A = 1.0
k = scc.Boltzmann
T = 10454
B = 1/(k*T)
L = 1*10**(-8)
eee = 10.22
sigma = 2.556*(10**(-10))
atoms1 = [[0,0,0]]
atoms2 = [[L,L,0]]
atoms3 = [[L,0,L]]
atoms4 = [[0,L,L]]
atoms1new =[0,0,0]
atoms2new =[L,L,0]
atoms3new =[L,0,L]
atoms4new =[0,L,L]
def roll_for_atom():
    a = rd.randint(1,4)
    if a == 1 :
        b = atoms1new
    elif (a==2) :
        b = atoms2new
    elif( a ==3):
        b = atoms3new
    else:
        b = atoms4new
    return b
def E(x,y,z):
    r1 = ((x-xi)**2+(y-yi)**2+(z-zi)**2)**(0.5)
    r2 = ((x-xj)**2+(y-yj)**2+(z-zj)**2)**(0.5)
    r3 = ((x-xk)**2+(y-yk)**2+(z-zk)**2)**(0.5)    
    Ee1 = 4*eee*(sigma**(12)/(r1**(12))-sigma**(6)/(r1**(6)))
    Ee2 = 4*eee*(sigma**(12)/(r2**(12))-sigma**(6)/(r2**(6)))
    Ee3 = 4*eee*(sigma**(12)/(r3**(12))-sigma**(6)/(r3**(6)))
    E = Ee1+Ee2+Ee3
    return E
j = 0
while j < 10000 :
    The_thing = roll_for_atom()        
    a = atoms1new
    b = atoms2new
    c = atoms3new
    d = atoms4new
    e = [a,b,c,d]
    e.remove(The_thing)
    xi =e[0][0] 
    yi =e[0][1]
    zi =e[0][2]
    xj =e[1][0]
    yj =e[1][1]
    zj =e[1][2]
    xk =e[2][0]
    yk =e[2][1]
    zk =e[2][2]
    
    x0 = The_thing[0]
    y0 = The_thing[1]
    z0 = The_thing[2]
    
    x = (rd.random())*L   #0~~L  
    y = (rd.random())*L
    z = (rd.random())*L    
    ei = E(x0,y0,z0)
    ef = E(x,y,z)
    de = ef - ei
    if de<=0 :
        if The_thing == atoms1new:
            atoms1.append([x,y,z])
            atoms1new = [x,y,z]
        elif The_thing == atoms2new:
            atoms2.append([x,y,z])
            atoms2new = [x,y,z]
        elif The_thing == atoms3new:
            atoms3.append([x,y,z])
            atoms3new = [x,y,z]
        elif The_thing == atoms4new:
            atoms4.append([x,y,z])
            atoms4new = [x,y,z]
 
    else :
        r = rd.random()
        re = math.exp(-(B)*de)
        if r <= re:
            if The_thing == atoms1new:
                atoms1.append([x,y,z])
                atoms1new = [x,y,z]
            elif The_thing == atoms2new:
                atoms2.append([x,y,z])
                atoms2new = [x,y,z]
            elif The_thing == atoms3new:
                atoms3.append([x,y,z])
                atoms3new = [x,y,z]
            elif The_thing == atoms4new:
                atoms4.append([x,y,z])
                atoms4new = [x,y,z]
            
        else: 
            if The_thing == atoms1new:
                atoms1.append([x0,y0,z0])
                atoms1new = [x0,y0,z0]
            elif The_thing == atoms2new:
                atoms2.append([x0,y0,z0])
                atoms2new = [x0,y0,z0]
            elif The_thing == atoms3new:
                atoms3.append([x0,y0,z0])
                atoms3new = [x0,y0,z0]
            elif The_thing == atoms4new:
                atoms4.append([x0,y0,z0])
                atoms4new = [x0,y0,z0]
    j+=1
fullofdots = atoms1+atoms2+atoms3+atoms4

j = 0
N = len(fullofdots)
xx = []
yy = []
zz = []    
print(N)
while j < N :
    xx.append(fullofdots[j][0])
    yy.append(fullofdots[j][1])
    zz.append(fullofdots[j][2])
    j+=1
print len(xx)
print len(yy)
print len(zz)



ax = plt.subplot(111,projection='3d')
#ax.hist((xx,yy,zz))
ax.plot(xx,yy,zz,'o')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()