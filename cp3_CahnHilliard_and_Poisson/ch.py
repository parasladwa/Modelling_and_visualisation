import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d













def initial_condition(N, case):
    
    arr = np.zeros((N, N))
    
    noise = np.random.rand(N, N)*0.1
    arr += noise
    
    if case == 'zeros':
        return arr
    
    elif case == 'halfs':
        arr = np.random.choice([-0.5, 0.5], size=(N, N))
        
    elif case == 'half':
        arr = np.random.choice([0.5, 0], size = (N, N))

    
    
    #randomly selected array of 0.5 or -0.5
    # arr = np.random.choice([0.5, -0.5], size=(N, N))
    # arr += noise
    return arr










#questions
#     free energy where it is? or before update phi

#to do:
#argparse initial conditions



def main():
        
    parser = argparse.ArgumentParser(description="Solve the Cahn-Hilliard equation")

    parser.add_argument("-a", "--a", type=float, default=0.1, help="Parameter a")
    parser.add_argument("-M", "--M", type=float, default=0.1, help="Parameter M")
    parser.add_argument("-k", "--kappa", type=float, default=0.1, help="Parameter kappa")
    
    parser.add_argument("-N", "--N", type=int, default=100, help="Grid size")
    parser.add_argument("-dt", "--dt", type=float, default=1, help="Time step")
    parser.add_argument("-dx", "--dx", type=float, default=1, help="Space step")
    
    parser.add_argument("-ns", "--num_steps", type=int, default=10000, help="Number of steps")
    parser.add_argument("-nth", "--show_nth", type=int, default=100, help="Show every nth step")
    parser.add_argument("-s", "--show_anim", action="store_true", help="show animation")
    
    parser.add_argument("-c", "--case", type=str, default = 'zeros', help = "choose a case from {zeros, halfs}")
    args = parser.parse_args()
    
    N = args.N
    a = args.a
    M = args.M
    kappa = args.kappa
    dt = args.dt
    dx = args.dx
    num_steps = args.num_steps
    show_nth = args.show_nth
    show_anim = args.show_anim
    
    
    
    
    
    
    kernel = np.array([[0, 1, 0],
                        [1, 0, 1],
                        [0, 1, 0]])
    
    kernel_j = np.array([[0, 0, 0],
                        [-1, 0, 1],
                        [0, 0, 0]])
    
    kernel_i = np.array([[0, -1, 0],
                        [0, 0, 0],
                        [0, 1, 0]])
    

    
    
    
    phi = initial_condition(N, args.case)
   
   
    free_energy_densities = np.zeros(num_steps)
   
   
   
    s = time.time()
    for step in range(num_steps): 
        
        
        #sim steps
        conv_phi= convolve2d(phi, kernel, mode='same', boundary='wrap')
        
        mu = -a*phi + a*phi**3 - kappa*(conv_phi - 4*phi)/(dx**2)
        
        conv_mu = convolve2d(mu, kernel, mode='same', boundary='wrap')
        phi_new = phi +  ((M*dt)/(dx**2)) * (conv_mu - 4*mu)
                

        conv_i = convolve2d(phi, kernel_i, mode='same', boundary='wrap')
        conv_j = convolve2d(phi, kernel_j, mode='same', boundary='wrap')
        
        
        t1 = -a * (phi**2)/2
        t2 = a * (phi**4)/4
        t3 = (kappa/2) * ((conv_i**2 + conv_j**2) /(2*dx)**2)

        free_energy = np.sum(t1 + t2 + t3)
        free_energy_densities[step] = free_energy
        
        
        phi = np.copy(phi_new)
        if step % show_nth == 0 and show_anim:
            plt.cla()
            plt.imshow(phi, animated=True)
            plt.text(0, -5, f"Step: {step}")
            plt.draw()
            plt.pause(0.01)

    en = time.time()
    print(f"Time: {en-s}")
    
    
    
    plt.figure()
    plt.plot(free_energy_densities)
    plt.title(f"Free energy density, case : {args.case}")
    plt.xlabel("Step")
    plt.ylabel("Free energy")
    plt.show()
    
main()

    
    

