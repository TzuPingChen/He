# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt


k=1.3806448E-23
def Cv(beta):
    return 3.0*k/2+1.734*k*beta**(-0.486)

xi = np.linspace(0.01,10,1000)
plt.subplot(221)
plt.plot(xi,Cv(xi))
plt.xlabel(r'$\beta$')
plt.ylabel("$C_V$("+r'$\beta)$') 
plt.subplot(222)
plt.plot(1.0/k/xi,Cv(xi))
plt.xlabel('T')
plt.ylabel("$C_V$(T)") 
plt.show()
