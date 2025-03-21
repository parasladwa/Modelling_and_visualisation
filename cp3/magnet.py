import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve





def main():
    
    dx = 0
    
    parser = argparse.ArgumentParser(description="Solve the Cahn-Hilliard equation")

    parser.add_argument("-N", "--N", type=int, default=100, help="Grid size")
    parser.add_argument("-nth", "--show_nth", type=int, default=100, help="Show every nth step")
    parser.add_argument("-s", "--show_anim", action="store_true", help="show animation")
    parser.add_argument("-t", "--termination_condition", type = float, default = 0.001, help = 'termination condition')

    
    args = parser.parse_args()
    
    j = np.zeros((args.N, args.N, args.N))
    A = np.copy(j)
    j[0, 0, :] = 1
    
        
    #3d kernel
    kernel = np.zeros((3, 3, 3))
    kernel[0, 1, 1] = 1 
    kernel[2, 1, 1] = 1
    kernel[1, 0, 1] = 1  
    kernel[1, 2, 1] = 1  
    kernel[1, 1, 0] = 1  
    kernel[1, 1, 2] = 1 
    
    
    
    #steps conv
    start = time.time()
    diff = np.inf
    while diff > args.termination_condition:
        
        conv = convolve(A, kernel, mode = 'constant')
        prev_A = np.copy(A)
        
        A = (1/6)*(conv + j)
        
        diff = np.linalg.norm(A - prev_A)


        if args.show_anim and args.show_nth:
            plt.cla()
            plt.imshow(A[0], animated=True)
            plt.draw()
            plt.pause(0.01)

    print(f"time = {time.time() - start}")

    plt.figure()
    plt.imshow(A[0])
    plt.colorbar()
    plt.title("resulting potential")
    plt.show()

    
    
    
    
    
    
    
    #gauss seidel
    
    print('GAUS SEIDEL METHOD')
    
    #checkerboard
    
    
    #split even and odd locs
    evens = np.empty((0, 3), dtype=int)  
    odds = np.empty((0, 3), dtype=int)  
    
    for i in range(args.N):
        for j in range(args.N):
            for k in range(args.N):
                
                if (i+j+k) % 2 == 0:
                    evens = np.vstack([evens, [i, j, k]])
                
                else:
                    odds = np.vstack([odds, [i, j, k]]) 
                
                
    print(evens)


main()