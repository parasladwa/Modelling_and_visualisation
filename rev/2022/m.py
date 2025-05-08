import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
from matplotlib.colors import ListedColormap






def compute_type_field(a, b, c):
    type_field = np.empty(shape = (len(a), len(a)))
    
    
    for i in range(len(a)):
        for j in range(len(a)):
            
            vals = [1-a[i,j]-b[i,j]-c[i,j], a[i, j], b[i, j], c[i, j]]
            
            type_field[i, j] = np.argmax(vals)

    return type_field





def main():
    
    nsteps = 5000
    shownth = 10
    N = 50
    D = 1
    q = 1
    p = 0.5
    dx = 1
    dt = 0.1

    
    cmap = ListedColormap(['gray', 'red', 'green', 'blue'])
    
    
    kernel = [[0, 1, 0],
              [1, 0, 1],
              [0, 1, 0]]
    
    
    # conv_phi= convolve2d(phi, kernel, mode='same', boundary='wrap')
    # mu = -a*phi + a*phi**3 - kappa*(conv_phi - 4*phi)/(dx**2)

    a = np.random.rand(N, N)/3
    b = np.random.rand(N, N)/3
    c = np.random.rand(N, N)/3
    
    for step in range(nsteps):
        
        
        
        #updates
        abc = 1-a-b-c
        
        conv_a = convolve2d(a, kernel, mode = 'same', boundary = 'wrap')
        a += (D*(conv_a-4*a)/(dx**2) + q*a*(abc) - p*a*c)*dt
        
        conv_b = convolve2d(b, kernel, mode = 'same', boundary = 'wrap')
        b += (D*(conv_b-4*b)/(dx**2) + q*b*(abc) - p*a*b)*dt
        
        conv_c = convolve2d(c, kernel, mode = 'same', boundary = 'wrap')
        c += (D*(conv_c-4*c)/(dx**2) + q*c*(abc) - p*b*c)*dt
        
        
        
        type_field = compute_type_field(a, b, c)
        
        if step%shownth == 0: 
            plt.cla()
            plt.imshow(type_field, animated=True, cmap = cmap)
            plt.text(0, -.1*N, f"Step: {step}")
            plt.draw()
            plt.pause(0.01)
    
    
    
    
main()