import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve





def main():
    
    dx = 1
    step = 0
    
    parser = argparse.ArgumentParser(description="Solve the Cahn-Hilliard equation")

    parser.add_argument("-N", "--N", type=int, default=50, help="Grid size")
    parser.add_argument("-nth", "--show_nth", type=int, default=100, help="Show every nth step")
    parser.add_argument("-s", "--show_anim", action="store_true", help="show animation")
    parser.add_argument("-t", "--termination_condition", type = float, default = 0.001, help = 'termination condition')
    parser.add_argument("-m", "--method", type = str, default = "gs", help = "choose method from gs or j, (gauss siedel or jacobi)")
    parser.add_argument("-w", "--omega", type = float, default = 1.0, help = "SOR factor, w > 1")
    parser.add_argument("-a", "--auto_gs", action="store_true", help="sim across w vals")
    parser.add_argument("-l", "--log", action="store_true", help="log datafiles")
    
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
            outfile_gs = open('gs_data.txt.', 'w')
            outfile_gs.write(f"<omega> <steps>\n")
            omegas = np.linspace(1, 2, )
                        
        
        print('GAUS SEIDEL METHOD')
          
        
        
        
        
        
        #naive approach

        
        
        for omega in omegas:
            step = 0
            start = time.time()
            print('omega', omega)
            Ap =  np.pad(A, pad_width=1, mode='constant', constant_values=0)
            diff = np.inf
            
            np.zeros((args.N, args.N, args.N))
            
            while diff > args.termination_condition and step < 2000:
                
                step+=1
                A_prev = np.copy(Ap)
                
                for i in range(1, args.N+1):
                    for jj in range(1, args.N+1):
                        for k in range(1, args.N+1):
                            
                            iterm = Ap[i+1, jj, k] + Ap[i-1, jj, k]
                            jterm = Ap[i, jj+1, k] + Ap[i, jj-1, k]
                            kterm = Ap[i, jj, k+1] + Ap[i, jj, k-1]
                            
                            ##over relaxation
                            new = (1/6) * (iterm+jterm+kterm + (dx**2)*j[i-1, jj-1, k-1])
                            
                            delta = new - Ap[i, jj, k]
                            
                            Ap[i, jj, k]= Ap[i, jj, k] + omega * delta
                            
                            



                            
                
                diff = np.linalg.norm(Ap[1:-1, 1:-1, 1:-1] - A_prev[1:-1, 1:-1, 1:-1])
                print(diff, step)

                
                

                        
                if args.show_anim:
                
                    plt.cla()
                    plt.imshow(Ap[1:-1, 1:-1, 1:-1][args.N//2], animated=True)
                    plt.title(f"step {step}")
                    plt.draw()
                    plt.pause(0.01)        

            
            if args.auto_gs:
                print('TIME', time.time()-start)
                outfile_gs.write(f"{omega} {step}\n")

            
            
        if not args.auto_gs:
            plt.figure()            
            plt.cla()
            plt.imshow(Ap[1:-1, 1:-1, 1:-1][args.N//2], animated=True)
            plt.show()
            
    # B field
    if not args.auto_gs:
        
        if args.method == 'gs':
            A = Ap[1:-1, 1:-1, 1:-1]
        
        
        grad_x_kernel = np.zeros((3, 3, 3))
        grad_y_kernel = np.zeros((3, 3, 3))
        grad_z_kernel = np.zeros((3, 3, 3))
        
        grad_x_kernel[2, 1, 1] = 1
        grad_x_kernel[0, 1, 1] = -1
        
        grad_y_kernel[1, 2, 1] = 1
        grad_y_kernel[1, 0, 1] = -1
        
        grad_z_kernel[1, 1, 2] = 1
        grad_z_kernel[1, 1, 0] = -1
        
        
        
        
        grad_x =  ( convolve(A, grad_x_kernel, mode='constant') ) / (2 * dx ) 
        grad_y =  ( convolve(A, grad_y_kernel, mode='constant'))  / (2 * dx ) 
        grad_z =  ( convolve(A, grad_z_kernel, mode='constant') ) / (2 * dx )
        
        Bx = np.copy(grad_y)
        By = -np.copy(grad_x)
        Bz = np.copy(grad_z)
        
        coords = np.arange(0, args.N)
        
        
        
        magnitudes = ( grad_x[:, :, args.N//2]**2 + grad_y[:, :, args.N//2]**2 ) ** (1/2)
        magnitudes *= 1.5
        
        plt.figure()
        plt.quiver(coords, coords, By[:, :, args.N//2]/magnitudes, Bx[:, :, args.N//2]/magnitudes)
        plt.title("B field")
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
            
            
            
    #outfiles
    if args.method == 'j' and args.log:
        outfile_j = open('magnet_jacobi.txt.', 'w')
        outfile_j.write(f"<i> <j> <vector_potential> <magnetic_field>\n")
        
        A_slice = A[args.N//2]
        
        Bx_s = Bx[:, :, args.N//2]
        By_s = By[:, :, args.N//2]
        Bz_s = Bz[:, :, args.N//2]
        
        for i in range(len(A[args.N//2])):
            for j in range(len(A[args.N//2, i])):
                outfile_j.write(f"{i} {j} {A_slice[i, j]} {Bx_s[i, j]} {By_s[i, j]} {Bz_s[i, j]}\n")

    
    if args.method == 'gs' and args.log:
        outfile_gs = open('magnet_gs', 'w')
        outfile_gs.write(f"<i> <j> <vector_potential> <magnetic_field>\n")
        
        A_slice = A[args.N//2]
        
        Bx_s = Bx[:, :, args.N//2]
        By_s = By[:, :, args.N//2]
        Bz_s = Bz[:, :, args.N//2]
        
        for i in range(len(A[args.N//2])):
            for j in range(len(A[args.N//2, i])):
                outfile_j.write(f"{i} {j} {A_slice[i, j]} {Bx_s[i, j]} {By_s[i, j]} {Bz_s[i, j]}\n")
        
    
    

main()






