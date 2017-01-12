import matplotlib.pyplot as plt  
import mpl_toolkits.mplot3d
import random
import numpy as np
import math
from scipy.spatial.distance import cdist
from visual import*
from visual import*
import sys
from types import*
from time import clock , time
import cPickle as pk

x=22
a=[12]
b='Bad apple'
x_=str(x)+'.txt'
print x_

f=file(x_,'w') 
pk.dump(a,f)
f.close()
f=file(x_)
print pk.load(f)

