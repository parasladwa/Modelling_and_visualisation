import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve





def main():
    
    dx = 1
    step = 0
    
    parser = argparse.ArgumentParser(description="Solve the Cahn-Hilliard equation")

    parser.add_argument("-N", "--N", type=int, default=100, help="Grid size")
    parser.add_argument("-nth", "--show_nth", type=int, default=100, help="Show every nth step")
    parser.add_argument("-s", "--show_anim", action="store_true", help="show animation")
    parser.add_argument("-t", "--termination_condition", type = float, default = 0.01, help = 'termination condition')
    parser.add_argument("-m", "--method", type = str, default = "gs", help = "choose method from gs or j, (gauss siedel or jacobi)")
    parser.add_argument("-w", "--omega", type = float, default = 1.0, help = "SOR factor, w > 1")
    parser.add_argument("-a", "--auto_gs", action="store_true", help="sim across w vals")

    args = parser.parse_args()
    
    j = np.zeros((args.N, args.N, args.N))
    A = np.copy(j)
    j[args.N//2, args.N//2, :] = 1
    
        
    #3d kernel
    kernel = np.zeros((3, 3, 3))
    kernel[0, 1, 1] = 1 
    kernel[2, 1, 1] = 1
    kernel[1, 0, 1] = 1  
    kernel[1, 2, 1] = 1  
    kernel[1, 1, 0] = 1  
    kernel[1, 1, 2] = 1 
    
    
    
    #steps conv
    if args.method == 'j':
        
        
        outfile_j = open('outfile_jacobi.txt', 'w')
        
        
        
        start = time.time()
        diff = np.inf
        while diff > args.termination_condition:
            
            step+=1
            
            conv = convolve(A, kernel, mode = 'constant')
            prev_A = np.copy(A)
            
            A = (1/6)*(conv + j)
            
            diff = np.linalg.norm(A - prev_A)


            if args.show_anim and args.show_nth:
                plt.cla()
                plt.imshow(A[args.N//2], animated=True)
                plt.draw()
                plt.title(f"step {step}")
                plt.pause(0.01)
            print(diff)

        print(f"time = {time.time() - start}")
        
        
        
        

        plt.figure()
        plt.imshow(A[args.N//2])
        plt.colorbar()
        plt.title("resulting potential")
        plt.show()

    

    
    
    
    
    
    #gauss seidel
    
    elif args.method == 'gs':
        
        
        omegas = np.array([args.omega], dtype = float)
        
        
        if args.auto_gs:
            outfile_gs = open('gs_data_magnet.txt.', 'w')
            omegas = np.linspace(1.45, 2, 10)
        
                        
        
        print('GAUS SEIDEL METHOD')
          
        
        #naive approach
        Ap =  np.pad(A, pad_width=1, mode='constant', constant_values=0)
        diff = np.inf

        
        
        for omega in omegas:
            steps = 0
            start = time.time()

            
            while diff > args.termination_condition:
                
                step+=1
                A_prev = np.copy(Ap)
                
                for i in range(1, args.N+1):
                    for jj in range(1, args.N+1):
                        for k in range(1, args.N+1):
                            
                            iterm = Ap[i+1, jj, k] + Ap[i-1, jj, k]
                            jterm = Ap[i, jj+1, k] + Ap[i, jj-1, k]
                            kterm = Ap[i, jj, k+1] + Ap[i, jj, k-1]
                            
                            Ap[i, jj, k] = (1/6) * (iterm+jterm+kterm + (dx**2)*j[i-1, jj-1, k-1])
                            


                            
                
                diff = np.linalg.norm(Ap[1:-1, 1:-1, 1:-1] - A_prev[1:-1, 1:-1, 1:-1])
                print(diff)
                

                
                

                        
                if args.show_anim:
                
                    plt.cla()
                    plt.imshow(Ap[1:-1, 1:-1, 1:-1][args.N//2], animated=True)
                    plt.title(f"step {step}")
                    plt.draw()
                    plt.pause(0.01)        

            print(omega, steps)
            print(time.time-start)
            
            
        if not args.auto_gs:
            plt.figure()            
            plt.cla()
            plt.imshow(Ap[1:-1, 1:-1, 1:-1][0], animated=True)
            plt.show()
            
        
        
        


main()




