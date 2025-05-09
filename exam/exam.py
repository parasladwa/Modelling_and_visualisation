from scipy.signal import convolve2d
import matplotlib.pyplot as plt
import numpy as np
import argparse





def simulate(N= 50, phi_0=0.1, animation=True, show_nth=10, max_steps=5000):
    #simulation parameters
    M = 0.1
    q_0 = 0.5
    a=0.1
    
    dx =1 
    dt = 0.01
    

    
    phi = np.full((N, N), phi_0, dtype=float)
    # + random noise
    noise = np.random.uniform(-0.1, 0.1, size=(N, N))
    phi += noise
    
    #kernel for convolution
    kernel = [[0, 1, 0],
                [1, -4, 1],
                [0, 1, 0]]
    
    
    step=0
    while step < max_steps:
        step += 1
        
        #simulation steps with central finite difference scheme using convolution
        lap_phi = convolve2d(phi, kernel, mode='same', boundary='wrap')
        lap_lap_phi = convolve2d(lap_phi, kernel, mode='same', boundary='wrap')
        
        mu = -a*phi + phi**3 + (q_0**4) * phi + 2*(q_0**2) * (lap_phi) + (lap_lap_phi)
        
        lap_mu = convolve2d(mu, kernel, mode='same', boundary='wrap')
        phi_new = phi + (M*(lap_mu-4*mu) *dt/(dx**2))
        
        phi = np.copy(phi_new)
        
        
        #animation
        if step % show_nth == 0 and animation:
            plt.cla()
            plt.title(f"Phase Field Simulation where phi_0 = {phi_0}")
            plt.imshow(phi, animated=True)
            plt.text(0, -.1*N, f"Step: {step}")
            plt.draw()
            plt.pause(0.01)
            print(np.mean(phi**2) - np.mean(phi)**2)
            
        #variance computation
        variance = np.mean(phi**2) - np.mean(phi)**2
        
        if step == max_steps:
            print(f"Final variance: {variance}, {step} steps, phi_0 = {phi_0}")
            return variance




        
def part_c():
    # Part c: Steady state variance vs initial phi_0
    # run simulation for different initial phi_0 values and plot the variance
    # for each initial phi_0, run the simulation and compute the variance
    # assume equlilibrium is reached after 20000 steps
    
    phi_0s = np.linspace(0.0, 0.25, 20)
    max_steps = 20000
    variances = np.zeros(len(phi_0s))
    
    for i in range(len(phi_0s)):
        phi_0 = phi_0s[i]
        variances[i] = simulate(phi_0=phi_0, animation=False, max_steps=max_steps)
        
    #open output file
    with open("part_c_output.txt", "w") as f:
        f.write("phi_0 variance\n")
        for i in range(len(phi_0s)):
            f.write(f"{phi_0s[i]} {variances[i]}\n")
    
    plt.figure()
    plt.scatter(phi_0s, variances)
    plt.xlabel("Initial phi_0")
    plt.ylabel("Variance")
    plt.title("steady state variance vs initial phi_0")
    plt.show()












def simulate_part_d(N= 50, phi_0=0.1, animation=True, show_nth=10, max_steps=1000, v_0=0.1):
    # Part d: Advection
    
    M = 0.1
    q_0 = 0.5
    a=0.1
    
    dx =1 
    dt = 0.01
        
    phi = np.full((N, N), phi_0, dtype=float)
    # + random noise
    noise = np.random.uniform(-0.1, 0.1, size=(N, N))
    phi += noise
    
    kernel = [[0, 1, 0],
                [1, -4, 1],
                [0, 1, 0]]
    
    
    # velocity field
    ys = np.arange(N)
    vx = -v_0 * np.sin(2 * np.pi * ys / N)  # shape (N,)
    
    # central difference kernel for x-derivative
    kernel_dx = np.array([[0, 0, 0],
                          [-0.5, 0, 0.5],
                          [0, 0, 0]])
    
    
    step=0
    while step < max_steps:
        step += 1
        
        #simulation steps with central finite difference scheme using convolution
        lap_phi = convolve2d(phi, kernel, mode='same', boundary='wrap')
        lap_lap_phi = convolve2d(lap_phi, kernel, mode='same', boundary='wrap')
        
        mu = -a*phi + phi**3 + (q_0**4) * phi + 2*(q_0**2) * (lap_phi) + (lap_lap_phi)
        
        lap_mu = convolve2d(mu, kernel, mode='same', boundary='wrap')

        # Advection term - only need x since vy = 0
        dphi_dx = convolve2d(phi, kernel_dx, mode='same', boundary='wrap')
        
        
        advect = np.zeros_like(phi)
        for i in range(N):
            for j in range(N):
                advect[i, j] = -vx[i] * dphi_dx[i, j]
        
        
        phi_new = phi + dt * (M * (lap_mu - 4 * mu)  + advect)
        phi = np.copy(phi_new)
    
        
        #animation
        if step % show_nth == 0 and animation:
            plt.cla()
            plt.title(f"Phase Field Simulation with advection {v_0} where phi_0 = {phi_0}")
            plt.imshow(phi, animated=True)
            plt.text(0, -.1*N, f"Step: {step}")
            plt.draw()
            plt.pause(0.01)
        
        variance = np.mean(phi**2) - np.mean(phi)**2
        
        if step == max_steps:
            print(f"Final variance: {variance}, {step} steps")
            return variance










def main():
    
    
    parser = argparse.ArgumentParser(description="Phase Field Simulation")
    parser.add_argument("--part", type=str, choices=["a", "b", "c", "d"], required=True, help="Part of the exam to run")
    parser.add_argument("--phi_0", type=float, default=0.1, help="Initial value of phi")
    parser.add_argument("--animation", action="store_true", help="Show animation")
    parser.add_argument("--show_nth", type=int, default=10, help="Show every nth step") 
    parser.add_argument("--max_steps", type=int, default=1000, help="Maximum number of steps")
    parser.add_argument("--N", type=int, default=50, help="Size of the grid")
    parser.add_argument("--v_0", type=float, default=0.1, help="Velocity for advection")
     
    args = parser.parse_args()
    
    if args.part == "a":
        simulate(N = args.N, phi_0=args.phi_0, animation=True, show_nth=args.show_nth, max_steps=args.max_steps)

    elif args.part == "b":
        print(f"refer to images stored in this directory for part b, and README.txt for discussion")
    
    elif args.part == "c":
        part_c()
    
    elif args.part == "d":
        simulate_part_d(N = args.N, animation=True, show_nth=args.show_nth, max_steps=args.max_steps, v_0=args.v_0)
        
       
       
       
       
main()






























