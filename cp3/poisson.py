import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve





def initial_condition(N):
    arr = np.zeros((N, N, N))
    return arr




def main():    
     
    dx = 1
    epsilon_zero = 1
    
    parser = argparse.ArgumentParser(description="Solve the Cahn-Hilliard equation")

    parser.add_argument("-N", "--N", type=int, default=100, help="Grid size")
    parser.add_argument("-nth", "--show_nth", type=int, default=100, help="Show every nth step")
    parser.add_argument("-s", "--show_anim", action="store_true", help="show animation")
    parser.add_argument("-t", "--termination_condition", type = float, default = 0.001, help = 'termination condition')
    
    args = parser.parse_args()
    
    
    
    phi = initial_condition(args.N)
    
    #point charge
    rho = initial_condition(args.N)
    rho[args.N//2, args.N//2, args.N//2] = 1

    
    #3d kernel
    kernel = np.zeros((3, 3, 3))
    kernel[0, 1, 1] = 1 
    kernel[2, 1, 1] = 1
    kernel[1, 0, 1] = 1  
    kernel[1, 2, 1] = 1  
    kernel[1, 1, 0] = 1  
    kernel[1, 1, 2] = 1 
    
    
    
    
    #steps
    start = time.time()
    diff = np.inf
    while diff > args.termination_condition:
        
        conv = convolve(phi, kernel, mode = 'constant')
        previous_phi = np.copy(phi)
        
        phi = (1/6)*(conv + rho)
        
        diff = np.linalg.norm(phi - previous_phi)
    
        if args.show_anim and args.show_nth:
            plt.cla()
            plt.imshow(phi[args.N//2], animated=True)
            plt.draw()
            plt.pause(0.01)

    print(f"time = {time.time() - start}")
    
    
    #resulting potential
    plt.figure()
    plt.imshow(phi[args.N//2])
    plt.colorbar()
    plt.title("resulting potential")
    plt.show()
    
    
    
    
    #resulting E field

    # Define 3D gradient kernels for each axis
    grad_x_kernel = np.zeros((3, 3, 3))
    grad_y_kernel = np.zeros((3, 3, 3))
    grad_z_kernel = np.zeros((3, 3, 3))
    
    grad_x_kernel[2, 1, 1] = 1
    grad_x_kernel[0, 1, 1] = -1
    
    grad_y_kernel[1, 2, 1] = 1
    grad_y_kernel[1, 0, 1] = -1
    
    grad_z_kernel[1, 1, 2] = 1
    grad_z_kernel[1, 1, 0] = -1


    grad_x =  ( convolve(phi, grad_x_kernel, mode='constant') ) / (2 * dx ) 
    grad_y =  ( convolve(phi, grad_y_kernel, mode='constant'))  / (2 * dx ) 
    grad_z =  ( convolve(phi, grad_z_kernel, mode='constant') ) / (2 * dx )
    
    print(grad_x.shape)
    
    coords = np.arange(0, args.N)
    
    magnitudes = ( grad_x[:, :, args.N//2]**2 + grad_y[:, :, args.N//2]**2 ) ** (1/2)
    magnitudes *= 1.5

    
    plt.figure()
    plt.quiver(coords, coords, grad_y[:, :, args.N//2]/magnitudes, grad_x[:, :, args.N//2]/magnitudes)
    plt.title('Electric Field (2D Slice)')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()
    
    
    
    
    #outfile
    phi_slice = phi[args.N//2]
    g_x_s = grad_x[:, :, args.N//2]
    g_y_s = grad_y[:, :, args.N//2]
    g_z_s = grad_z[:, :, args.N//2]

    
    name = "poisson.txt"
    outfile = open(name, 'w')
    outfile.write(f"i j potential Ex Ey Ez\n")
    
    
    for i in range(len(phi[args.N//2])):
        for j in range(len(phi[args.N//2, i])):
            outfile.write(f"{i} {j} {phi_slice[i, j]} {g_x_s[i, j]}, {g_y_s[i, j]} {g_z_s[i, j]}\n")
            
            
    outfile.close()


main()
    
