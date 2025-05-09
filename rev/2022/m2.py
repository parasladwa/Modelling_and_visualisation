import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import time
from scipy.signal import convolve2d




def compute_type_field(a, b, c):
    
    tf = np.zeros((len(a), len(a[0])))
    
    for i in range(len(a)):
        for j in range(len(a[0])):
            
            l = [None, a[i, j], b[i, j], c[i, j]]
            l[0] = 1 - l[1] - l[2] -l[3]
            
            tf[i, j] = np.argmax(l)
    return tf




def simulate():
    
    
    kernel = [[0, 1, 0],
            [1, 0, 1],
            [0, 1, 0]]


    N = 50
    dx = 1
    
    dt = 0.1
    D = 1
    q = 1
    p = 0.5
    
    converged = False
        
    a = np.random.rand(N, N)/3
    b = np.random.rand(N, N)/3
    c = np.random.rand(N, N)/3

    nstep = 0
    
    while not converged:
        nstep+=1
        #steps
        a_conv = convolve2d(a, kernel, mode = 'same', boundary = 'wrap')
        a_new = a + (D * (a_conv -4*a)+ q*a*(1-a-b-c) - p*a*c)*dt
        
        b_conv = convolve2d(b, kernel, mode = 'same', boundary = 'wrap')
        b_new = b + (D*(b_conv-4*b) + q*b*(1-a-b-c) - p*a*b)*dt
        
        c_conv = convolve2d(c, kernel, mode = 'same', boundary = 'wrap')
        c_new = c + (D*(c_conv-4*c)+ q*c*(1-a-b-c) - p*b*c)*dt
        
        a = np.copy(a_new)
        b = np.copy(b_new)
        c = np.copy(c_new)



        type_field = compute_type_field(a, b, c)
        #check convergence  
        


        if nstep%50 == 0:
            plt.cla()
            plt.imshow(type_field, animated=True)
            plt.draw()
            plt.pause(0.01)
        

    
    
    

simulate()





























