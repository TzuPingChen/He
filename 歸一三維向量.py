import matplotlib.pyplot as plt  
import mpl_toolkits.mplot3d
import random
import numpy as np
import math
from scipy.spatial.distance import cdist

def vec():
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
    return x,y,z

print vec()
