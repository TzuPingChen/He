# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

e=10.22
sig=0.2556*10**(-9)
Vo = 0.5*sig**3
sig6 = sig**6

def Utot(V):
    return 3*e*sig6*(((1.0/8)*sig6*V**(-4))-V**(-2))

def Pv(V):
    return 3*e*sig6*(((0.5)*sig6*V**(-5))-2*V**(-3))

def KT(V):
    return 1.0/(18*e*(sig6)*((5.0/12)*(sig6)*V**(-5)-V**(-3)))

xi = np.linspace(0.5*Vo,Vo,1000)
plt.subplot(223)
plt.plot(xi/Vo,KT(xi))
plt.xlabel("V(1/Vo)")
plt.ylabel("$\kappa_T(V)$") 
plt.subplot(222)
plt.plot(xi/Vo,Pv(xi))
plt.xlabel("V(1/Vo)")
plt.ylabel("P(V)")
plt.subplot(221)
plt.plot(xi/Vo,Utot(xi))
plt.xlabel("V(1/Vo)")
plt.ylabel("Utot(V)")
plt.show()
