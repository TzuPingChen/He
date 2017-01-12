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


a=[12]
b='Bad apple'
f=file('Nhe4,V.txt','w') 
pk.dump(a,f)
f.close()
f=file('Nhe4,V.txt')
print pk.load(f)
