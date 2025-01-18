import sys
import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import time



def set_up(N):
    #randomly set up array
    arr = np.zeros([N, N], dtype=int)
    
    for i in range(0, len(arr)):
        for j in range(0, len(arr[0])):
            arr[i][j] = random.choice([-1, 1])
    
    return arr


arr = set_up(50)
#copied
fig = plt.figure()
im=plt.imshow(arr, animated=True)

for i in range(1000):
    arr = set_up(50)
  
    plt.cla()
    im=plt.imshow(arr, animated=True)
    plt.draw()
    plt.pause(0.0001)

    
    
    
    
    
    
    
    
    