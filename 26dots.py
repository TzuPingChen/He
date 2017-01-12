# -*- coding: utf-8 -*-
import numpy as np
import random as rd
import matplotlib.pyplot as plt
import scipy.constants as scc
import math
import mpl_toolkits.mplot3d

a = 2
A = 1.0
k = scc.Boltzmann
T = 200
B = 1/(k*T)
L = 5*10**(-6)
eee = 10.22
sigma = 2.556*(10**(-10))
atoms1 = [[0,0,0]]
atoms2 = [[L,L,0]]
atoms3 = [[L,0,L]]
atoms4 = [[0,L,L]]
atoms5 = [[L,L,L]]
atoms6 = [[0,0,L]]
atoms7 = [[0,L,0]]
atoms8 = [[L,0,0]]
atoms9 = [[L/2,L,L]]
atoms10 = [[L,L/2,L]]
atoms11 = [[L,L,L/2]]
atoms12 = [[L/2,L/2,L]]
atoms13 = [[L,L/2,L/2]]
atoms14 = [[L/2,L,L/2]] 
atoms15 = [[0,L,L/2]]
atoms16 = [[0,L/2,0]]
atoms17 = [[L/2,0,0]]
atoms18 = [[0,0,L/2]]
atoms19 = [[0,L/2,L/2]]
atoms20 = [[L/2,L/2,0]]
atoms21 = [[L/2,0,L/2]]
atoms22 = [[0,L/2,L]]
atoms23 = [[L/2,L,0]]
atoms24 = [[L/2,0,L]]
atoms25 = [[L,L/2,0]]
atoms26 = [[L,0,L/2]]

atoms1new =[0,0,0]
atoms2new =[L,L,0]
atoms3new =[L,0,L]
atoms4new =[0,L,L]
atoms5new = [L,L,L]
atoms6new = [0,0,L]
atoms7new = [0,L,0]
atoms8new = [L,0,0]
atoms9new = [L/2,L,L]
atoms10new = [L,L/2,L]
atoms11new = [L,L,L/2]
atoms12new = [L/2,L/2,L]
atoms13new = [L,L/2,L/2]
atoms14new = [L/2,L,L/2] #
atoms15new = [L/2,L/2,L/2]
atoms16new = [0,L/2,0]
atoms17new = [L/2,0,0]
atoms18new = [0,0,L/2]
atoms19new = [0,L/2,L/2]
atoms20new = [L/2,L/2,0]
atoms21new = [L/2,0,L/2]
atoms22new = [0,L/2,L]
atoms23new = [L/2,L,0]
atoms24new = [L/2,0,L]
atoms25new = [L,L/2,0]
atoms26new = [L,0,L/2]
whole = atoms1+atoms2+atoms3+atoms4+atoms5+atoms6+atoms7+atoms8+atoms9+atoms10+atoms11+atoms12+atoms13+atoms14+atoms15+atoms16+atoms17+atoms18+atoms19+atoms20+atoms21+atoms22+atoms23+atoms24+atoms25+atoms26
xiii=[]
yiii=[]
ziii=[]
NN = len(whole)
cc=0
while cc<NN:
    xiii.append(whole[cc][0])
    yiii.append(whole[cc][1])
    ziii.append(whole[cc][2])
    cc+=1
def roll_for_atom():     #挑任一粒子位置輸出
    a = rd.randint(1,26)
    if a == 1 :
        b = atoms1new
    elif (a==2) :
        b = atoms2new
    elif( a ==3):
        b = atoms3new
    elif( a ==4):
        b = atoms4new
    elif( a ==5):
        b = atoms5new
    elif( a ==6):
        b = atoms6new
    elif( a ==7):
        b = atoms7new
    elif( a ==8):
        b = atoms8new
    elif a == 9 :
        b = atoms9new
    elif (a==10) :
        b = atoms10new
    elif( a ==11):
        b = atoms11new
    elif( a ==12):
        b = atoms12new
    elif( a ==13):
        b = atoms13new
    elif( a ==14):
        b = atoms14new
    elif( a ==15):
        b = atoms15new
    elif( a ==16):
        b = atoms16new
    elif a == 17:                 
        b = atoms17new
    elif (a==18) :
        b = atoms18new
    elif( a ==19):
        b = atoms19new
    elif( a ==20):
        b = atoms20new
    elif( a ==21):
        b = atoms21new
    elif( a ==22):
        b = atoms22new
    elif( a ==23):
        b = atoms23new
    elif( a ==24):
        b = atoms24new
    elif( a ==25):
        b = atoms25new
    else:
        b = atoms26new
    return b

def E(x,y,z):
    r1 = ((x-xa)**2+(y-ya)**2+(z-za)**2)**(0.5)
    r2 = ((x-xb)**2+(y-yb)**2+(z-zb)**2)**(0.5)
    r3 = ((x-xc)**2+(y-yc)**2+(z-zc)**2)**(0.5)
    r4 = ((x-xd)**2+(y-yd)**2+(z-zd)**2)**(0.5)
    r5 = ((x-xe)**2+(y-ye)**2+(z-ze)**2)**(0.5)
    r6 = ((x-xf)**2+(y-yf)**2+(z-zf)**2)**(0.5)
    r7 = ((x-xg)**2+(y-yg)**2+(z-zg)**2)**(0.5)
    r8 = ((x-xh)**2+(y-yh)**2+(z-zh)**2)**(0.5)
    r9 = ((x-xi)**2+(y-yi)**2+(z-zi)**2)**(0.5)
    r10 = ((x-xj)**2+(y-yj)**2+(z-zj)**2)**(0.5)
    r11 = ((x-xk)**2+(y-yk)**2+(z-zk)**2)**(0.5)
    r12 = ((x-xl)**2+(y-yl)**2+(z-zl)**2)**(0.5)
    r13 = ((x-xm)**2+(y-ym)**2+(z-zm)**2)**(0.5)
    r14 = ((x-xn)**2+(y-yn)**2+(z-zn)**2)**(0.5)
    r15 = ((x-xo)**2+(y-yo)**2+(z-zo)**2)**(0.5)
    r16 = ((x-xp)**2+(y-yp)**2+(z-zp)**2)**(0.5)
    r17 = ((x-xq)**2+(y-yq)**2+(z-zq)**2)**(0.5)
    r18 = ((x-xr)**2+(y-yr)**2+(z-zr)**2)**(0.5)
    r19 = ((x-xs)**2+(y-ys)**2+(z-zs)**2)**(0.5)
    r20 = ((x-xt)**2+(y-yt)**2+(z-zt)**2)**(0.5)
    r21 = ((x-xu)**2+(y-yu)**2+(z-zu)**2)**(0.5)
    r22 = ((x-xv)**2+(y-yv)**2+(z-zv)**2)**(0.5)
    r23 = ((x-xw)**2+(y-yw)**2+(z-zw)**2)**(0.5)
    r24 = ((x-xx)**2+(y-yx)**2+(z-zx)**2)**(0.5)
    r25 = ((x-xy)**2+(y-yy)**2+(z-zy)**2)**(0.5)
    Ee1 = 4*eee*(sigma**(12)/(r1**(12))-sigma**(6)/(r1**(6)))
    Ee2 = 4*eee*(sigma**(12)/(r2**(12))-sigma**(6)/(r2**(6)))
    Ee3 = 4*eee*(sigma**(12)/(r3**(12))-sigma**(6)/(r3**(6)))
    Ee4 = 4*eee*(sigma**(12)/(r4**(12))-sigma**(6)/(r4**(6)))
    Ee5 = 4*eee*(sigma**(12)/(r5**(12))-sigma**(6)/(r5**(6)))
    Ee6 = 4*eee*(sigma**(12)/(r6**(12))-sigma**(6)/(r6**(6)))
    Ee7 = 4*eee*(sigma**(12)/(r7**(12))-sigma**(6)/(r7**(6)))
    Ee8 = 4*eee*(sigma**(12)/(r8**(12))-sigma**(6)/(r8**(6)))
    Ee9 = 4*eee*(sigma**(12)/(r9**(12))-sigma**(6)/(r9**(6)))
    Ee10 = 4*eee*(sigma**(12)/(r10**(12))-sigma**(6)/(r10**(6)))
    Ee11 = 4*eee*(sigma**(12)/(r11**(12))-sigma**(6)/(r11**(6)))                #
    Ee12 = 4*eee*(sigma**(12)/(r12**(12))-sigma**(6)/(r12**(6)))
    Ee13 = 4*eee*(sigma**(12)/(r13**(12))-sigma**(6)/(r13**(6)))
    Ee14 = 4*eee*(sigma**(12)/(r14**(12))-sigma**(6)/(r14**(6)))
    Ee15 = 4*eee*(sigma**(12)/(r15**(12))-sigma**(6)/(r15**(6)))
    Ee16 = 4*eee*(sigma**(12)/(r16**(12))-sigma**(6)/(r16**(6)))
    Ee17 = 4*eee*(sigma**(12)/(r17**(12))-sigma**(6)/(r17**(6)))
    Ee18 = 4*eee*(sigma**(12)/(r18**(12))-sigma**(6)/(r18**(6)))
    Ee19 = 4*eee*(sigma**(12)/(r19**(12))-sigma**(6)/(r19**(6)))
    Ee20 = 4*eee*(sigma**(12)/(r20**(12))-sigma**(6)/(r20**(6)))
    Ee21 = 4*eee*(sigma**(12)/(r21**(12))-sigma**(6)/(r21**(6)))
    Ee22 = 4*eee*(sigma**(12)/(r22**(12))-sigma**(6)/(r22**(6)))
    Ee23 = 4*eee*(sigma**(12)/(r23**(12))-sigma**(6)/(r23**(6)))
    Ee24 = 4*eee*(sigma**(12)/(r24**(12))-sigma**(6)/(r24**(6)))
    Ee25 = 4*eee*(sigma**(12)/(r25**(12))-sigma**(6)/(r25**(6)))
    E = Ee1+Ee2+Ee3+Ee4+Ee5+Ee6+Ee7+Ee8+Ee9+Ee10+Ee11+Ee12+Ee13+Ee14+Ee15+Ee16+Ee17+Ee18+Ee19+Ee20+Ee21+Ee22+Ee23+Ee24+Ee25
    return E

cc = 0
while cc < 75000 :
    The_thing = roll_for_atom()        
    a = atoms1new
    b = atoms2new
    c = atoms3new
    d = atoms4new
    e = atoms5new
    f = atoms6new
    g = atoms7new
    h = atoms8new
    i = atoms9new
    j = atoms10new
    k = atoms11new
    l = atoms12new
    m = atoms13new
    n = atoms14new
    o = atoms15new
    p = atoms16new
    q = atoms17new
    r = atoms18new
    s = atoms19new
    t = atoms20new
    u = atoms21new
    v = atoms22new
    w = atoms23new
    x = atoms24new
    y = atoms25new
    z = atoms26new
    
    lala = [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z]
    lala.remove(The_thing)
    
    
    xa =lala[0][0] 
    ya =lala[0][1]
    za =lala[0][2]
    xb =lala[1][0]
    yb =lala[1][1]
    zb =lala[1][2]
    xc =lala[2][0]
    yc =lala[2][1]
    zc =lala[2][2]
    xd =lala[3][0]
    yd =lala[3][1]
    zd =lala[3][2]
    xe =lala[4][0]
    ye =lala[4][1]
    ze =lala[4][2]
    xf =lala[5][0]
    yf =lala[5][1]
    zf =lala[5][2]
    xg =lala[6][0]
    yg =lala[6][1]
    zg =lala[6][2]
    xh =lala[7][0]
    yh =lala[7][1]
    zh =lala[7][2]
    xi =lala[8][0] 
    yi =lala[8][1]
    zi =lala[8][2]
    xj =lala[9][0]
    yj =lala[9][1]
    zj =lala[9][2]
    xk =lala[10][0]
    yk =lala[10][1]
    zk =lala[10][2]
    xl =lala[11][0]
    yl =lala[11][1]
    zl =lala[11][2]
    xm =lala[12][0]
    ym =lala[12][1]
    zm =lala[12][2]
    xn =lala[13][0]
    yn =lala[13][1]
    zn =lala[13][2]
    xo =lala[14][0]
    yo =lala[14][1]
    zo =lala[14][2]
    xp =lala[15][0]
    yp =lala[15][1]
    zp =lala[15][2]
    xq =lala[16][0] 
    yq =lala[16][1]
    zq =lala[16][2]
    xr =lala[17][0]
    yr =lala[17][1]
    zr =lala[17][2]
    xs =lala[18][0]
    ys =lala[18][1]
    zs =lala[18][2]
    xt =lala[19][0]
    yt =lala[19][1]
    zt =lala[19][2]
    xu =lala[20][0]
    yu =lala[20][1]
    zu =lala[20][2]
    xv =lala[21][0]
    yv =lala[21][1]
    zv =lala[21][2]
    xw =lala[22][0]
    yw =lala[22][1]
    zw =lala[22][2]
    xx =lala[23][0]
    yx =lala[23][1]
    zx =lala[23][2]
    xy =lala[24][0]
    yy =lala[24][1]
    zy =lala[24][2]

    
    
    x0 = The_thing[0]     #initial spot
    y0 = The_thing[1]
    z0 = The_thing[2]
    
    x = (rd.random())*L*2   #0~~L  new spot
    y = (rd.random())*L*2
    z = (rd.random())*L*2    
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
        elif The_thing == atoms5new:
            atoms5.append([x,y,z])
            atoms5new = [x,y,z]
        elif The_thing == atoms6new:
            atoms6.append([x,y,z])
            atoms6new = [x,y,z]
        elif The_thing == atoms7new:
            atoms7.append([x,y,z])
            atoms7new = [x,y,z]
        elif The_thing == atoms8new:
            atoms8.append([x,y,z])
            atoms8new = [x,y,z]
        elif The_thing == atoms9new:
            atoms9.append([x,y,z])
            atoms9new = [x,y,z]
        elif The_thing == atoms10new:
            atoms10.append([x,y,z])
            atoms10new = [x,y,z]
        elif The_thing == atoms11new:
            atoms11.append([x,y,z])
            atoms11new = [x,y,z]
        elif The_thing == atoms12new:
            atoms12.append([x,y,z])
            atoms12new = [x,y,z]
        elif The_thing == atoms13new:
            atoms13.append([x,y,z])
            atoms13new = [x,y,z]
        elif The_thing == atoms14new:
            atoms14.append([x,y,z])
            atoms14new = [x,y,z]
        elif The_thing == atoms15new:
            atoms15.append([x,y,z])
            atoms15new = [x,y,z]
        elif The_thing == atoms16new:
            atoms16.append([x,y,z])
            atoms16new = [x,y,z]
        elif The_thing == atoms17new:
            atoms17.append([x,y,z])
            atoms17new = [x,y,z]
        elif The_thing == atoms18new:
            atoms18.append([x,y,z])
            atoms18new = [x,y,z]
        elif The_thing == atoms19new:
            atoms19.append([x,y,z])
            atoms19new = [x,y,z]
        elif The_thing == atoms20new:
            atoms20.append([x,y,z])
            atoms20new = [x,y,z]
        elif The_thing == atoms21new:
            atoms21.append([x,y,z])
            atoms21new = [x,y,z]
        elif The_thing == atoms22new:
            atoms22.append([x,y,z])
            atoms22new = [x,y,z]
        elif The_thing == atoms23new:
            atoms23.append([x,y,z])
            atoms23new = [x,y,z]
        elif The_thing == atoms24new:
            atoms24.append([x,y,z])
            atoms24new = [x,y,z]
        elif The_thing == atoms25new:
            atoms25.append([x,y,z])
            atoms25new = [x,y,z]
        elif The_thing == atoms26new:
            atoms26.append([x,y,z])
            atoms26new = [x,y,z]
 
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
            elif The_thing == atoms5new:
                atoms5.append([x,y,z])
                atoms5new = [x,y,z]
            elif The_thing == atoms6new:
                atoms6.append([x,y,z])
                atoms6new = [x,y,z]
            elif The_thing == atoms7new:
                atoms7.append([x,y,z])
                atoms7new = [x,y,z]
            elif The_thing == atoms8new:
                atoms8.append([x,y,z])
                atoms8new = [x,y,z]
            elif The_thing == atoms9new:
                atoms9.append([x,y,z])
                atoms9new = [x,y,z]
            elif The_thing == atoms10new:
                atoms10.append([x,y,z])
                atoms10new = [x,y,z]
            elif The_thing == atoms11new:
                atoms11.append([x,y,z])
                atoms11new = [x,y,z]
            elif The_thing == atoms12new:
                atoms12.append([x,y,z])
                atoms12new = [x,y,z]
            elif The_thing == atoms13new:
                atoms13.append([x,y,z])
                atoms13new = [x,y,z]
            elif The_thing == atoms14new:
                atoms14.append([x,y,z])
                atoms14new = [x,y,z]
            elif The_thing == atoms15new:
                atoms15.append([x,y,z])
                atoms15new = [x,y,z]
            elif The_thing == atoms16new:
                atoms16.append([x,y,z])
                atoms16new = [x,y,z]
            elif The_thing == atoms17new:
                atoms17.append([x,y,z])
                atoms17new = [x,y,z]
            elif The_thing == atoms18new:
                atoms18.append([x,y,z])
                atoms18new = [x,y,z]
            elif The_thing == atoms19new:
                atoms19.append([x,y,z])
                atoms19new = [x,y,z]
            elif The_thing == atoms20new:
                atoms20.append([x,y,z])
                atoms20new = [x,y,z]
            elif The_thing == atoms21new:
                atoms21.append([x,y,z])
                atoms21new = [x,y,z]
            elif The_thing == atoms22new:
                atoms22.append([x,y,z])
                atoms22new = [x,y,z]
            elif The_thing == atoms23new:
                atoms23.append([x,y,z])
                atoms23new = [x,y,z]
            elif The_thing == atoms24new:
                atoms24.append([x,y,z])
                atoms24new = [x,y,z]
            elif The_thing == atoms25new:
                atoms25.append([x,y,z])
                atoms25new = [x,y,z]
            elif The_thing == atoms26new:
                atoms26.append([x,y,z])
                atoms26new = [x,y,z]
            
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
            elif The_thing == atoms5new:
                atoms5.append([x0,y0,z0])
                atoms5new = [x0,y0,z0]
            elif The_thing == atoms6new:
                atoms6.append([x0,y0,z0])
                atoms6new = [x0,y0,z0]
            elif The_thing == atoms7new:
                atoms7.append([x0,y0,z0])
                atoms7new = [x0,y0,z0]
            elif The_thing == atoms8new:
                atoms8.append([x0,y0,z0])
                atoms8new = [x0,y0,z0]
            elif The_thing == atoms9new:
                atoms9.append([x0,y0,z0])
                atoms9new = [x0,y0,z0]
            elif The_thing == atoms10new:
                atoms10.append([x0,y0,z0])
                atoms10new = [x0,y0,z0]
            elif The_thing == atoms11new:
                atoms11.append([x0,y0,z0])
                atoms11new = [x0,y0,z0]
            elif The_thing == atoms12new:
                atoms12.append([x0,y0,z0])
                atoms12new = [x0,y0,z0]
            elif The_thing == atoms13new:
                atoms13.append([x0,y0,z0])
                atoms13new = [x0,y0,z0]
            elif The_thing == atoms14new:
                atoms14.append([x0,y0,z0])
                atoms14new = [x0,y0,z0]
            elif The_thing == atoms15new:
                atoms15.append([x0,y0,z0])
                atoms15new = [x0,y0,z0]
            elif The_thing == atoms16new:
                atoms16.append([x0,y0,z0])
                atoms16new = [x0,y0,z0]
            elif The_thing == atoms17new:
                atoms17.append([x0,y0,z0])
                atoms17new = [x0,y0,z0]
            elif The_thing == atoms18new:
                atoms18.append([x0,y0,z0])
                atoms18new = [x0,y0,z0]
            elif The_thing == atoms19new:
                atoms19.append([x0,y0,z0])
                atoms19new = [x0,y0,z0]
            elif The_thing == atoms20new:
                atoms20.append([x0,y0,z0])
                atoms20new = [x0,y0,z0]
            elif The_thing == atoms21new:
                atoms21.append([x0,y0,z0])
                atoms21new = [x0,y0,z0]
            elif The_thing == atoms22new:
                atoms22.append([x0,y0,z0])
                atoms22new = [x0,y0,z0]
            elif The_thing == atoms23new:
                atoms23.append([x0,y0,z0])
                atoms23new = [x0,y0,z0]
            elif The_thing == atoms24new:
                atoms24.append([x0,y0,z0])
                atoms24new = [x0,y0,z0]
            elif The_thing == atoms25new:
                atoms25.append([x0,y0,z0])
                atoms25new = [x0,y0,z0]
            elif The_thing == atoms26new:
                atoms26.append([x0,y0,z0])
                atoms26new = [x0,y0,z0]
    cc+=1
fullofdots = atoms1+atoms2+atoms3+atoms4+atoms5+atoms6+atoms7+atoms8+atoms9+atoms10+atoms11+atoms12+atoms13+atoms14+atoms15+atoms16+atoms17+atoms18+atoms19+atoms20+atoms21+atoms22+atoms23+atoms24+atoms25+atoms26

cc = 0
N = len(fullofdots)
xxx = []
yyy = []
zzz = []    

while cc < N :
    xxx.append(fullofdots[cc][0])
    yyy.append(fullofdots[cc][1])
    zzz.append(fullofdots[cc][2])
    cc+=1
print len(xxx)
print len(yyy)
print len(zzz)

xfff=[]
yfff=[]
zfff=[]
fullofdotspart2= atoms1new+atoms2new+atoms3new+atoms4new+atoms5new+atoms6new+atoms7new+atoms8new+atoms9new+atoms10new+atoms11new+atoms12new+atoms13new+atoms14new+atoms15new+atoms16new+atoms17new+atoms18new+atoms19new+atoms20new+atoms21new+atoms22new+atoms23new+atoms24new+atoms25new+atoms26new
NN = len(fullofdotspart2)
print NN

def drange(start, stop, step):
    
    r = start
    while r < stop:
    	yield r
    	r += step
     
    return

for xyz in drange(0,NN,3) :
    xfff.append(fullofdotspart2[xyz])
    yfff.append(fullofdotspart2[xyz+1])
    zfff.append(fullofdotspart2[xyz+2])


ax = plt.subplot(111,projection='3d')
#ax.plot(xxx,yyy,zs=0,zdir='z',c='yellow')
#ax.plot(xxx,yyy,zzz,'o',c='yellow')
ax.plot(xiii,yiii,ziii,'o',c='red')
ax.plot(xfff,yfff,zfff,'o',c='blue')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
