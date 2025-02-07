import numpy as np
from scipy.signal import convolve2d
import matplotlib.pyplot as plt


def initialise(N):
    arr = np.random.randint(0, 2, (N, N))
    return arr
   
    



def main():
    
    #constants/parameters
    N = 50
    show_anim = True
    nsteps = 1000
    
    arr = initialise(N)
    
    
    if show_anim:
        fig = plt.figure()
        im=plt.imshow(arr, animated=True)
    


    
    
    kernel = np.array([[1, 1, 1],
                       [1, 0, 1],
                       [1, 1, 1]])
    

    
    for n in range(nsteps):
        
        convolution = convolve2d(arr, kernel, mode="same", boundary="wrap")

        for i in range(N):
            for j in range(N):
                
                
                # alive and neighbours < 2 : dead
                if arr[i, j] == 1 and convolution[i, j] < 2:
                    arr[i, j] = 0
                
                # alive and neighbours == 2, 3 : remains
                elif arr[i, j] == 1 and 2<= convolution[i, j] <= 3:
                    pass
                
                # alive and neighbours > 3 : dead
                elif arr[i, j] == 1 and convolution[i, j] > 3:
                    arr[i, j] = 0
                    
                elif arr[i, j] == 0 and convolution[i, j] == 3:
                    arr[i, j] = 1
                    
                
                else:
                    pass
                    #print('unassigned case')
                    
        if show_anim:           
            plt.cla()
            im=plt.imshow(arr, animated=True)
            plt.draw()
            plt.pause(0.001)
    
    

    




main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    