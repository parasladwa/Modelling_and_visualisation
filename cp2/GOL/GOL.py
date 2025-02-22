import sys
import time
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

start = time.time()



def initialise(N, method):
    
    if method == 'random':
        arr = np.random.randint(0, 2, (N, N), dtype=int)
    
    elif method == 'oscillator':
        arr = np.zeros((N, N), dtype=int)
        arr[N//2, N//2-1:N//2+2] = 1
    
    elif method == 'glider':
        arr = np.zeros((N, N), dtype=int)
        arr[N//2-1, N//2] = 1
        arr[N//2, N//2+1] = 1
        arr[N//2+1, N//2-1:N//2+2] = 1

    return arr



def main():
    
    print(f"<run> <GOL.py> <N> <method>"
          f"methods = random, oscillator, glider")
    
    #constants/parameters
    method = 'random'
    N = 50
    n_simulations = 1000
    consider_equilibrated = 10
    max_steps = 3000
    show_anim = False
    show_plot = True
    steps_to_equil = np.array([], dtype=int)
    filename = 'GOL_data.txt'
    nsteps_other = 200
    coms = np.array([], dtype=float)
    vels = np.array([], dtype=float)


    
    if len(sys.argv) > 1:
        N = int(sys.argv[1])
        method = sys.argv[2]
        n_simulations = 1
        show_anim = True
        show_plot = False
        filename = 'user.txt'
        
        if method not in ['random', 'oscillator', 'glider']:
            print("invalid method")
            return
        
        
        if method != 'random':
            consider_equilibrated = nsteps_other
    
    
    outfile = open(filename, 'w')
    outfile.write("run_no steps_to_equilibrate")
    
    kernel = np.array([[1, 1, 1],
                       [1, 0, 1],
                       [1, 1, 1]])
    
    arr = np.empty([N, N])
    
    
    print(f"running {n_simulations} simulation(s)"
          f" with {N}x{N} grid and {method} initialisation")
    
    for n_sim in range(n_simulations):
        
        evol = np.array([], dtype=int)
        arr = initialise(N, method)
        evol = np.append(evol, [np.sum(arr)])
        equilibrated = False
        step = 0
        
        while not equilibrated:
            
            conv = convolve2d(arr, kernel, mode="same", boundary="wrap")
            
            #GOL conditions
            for i in range(N):
                for j in range(N):
                    if arr[i, j] == 1 and conv[i, j] < 2:
                        arr[i, j] = 0
                    elif arr[i, j] == 1 and 2 <= conv[i, j] <= 3:
                        pass
                    elif arr[i, j] == 1 and conv[i, j] > 3:
                        arr[i, j] = 0
                    elif arr[i, j] == 0 and conv[i, j] == 3:
                        arr[i, j] = 1
                    else:
                        pass

            step += 1
            evol = np.append(evol, [np.sum(arr)])
            
            #if equilibrated
            if np.all(evol[-consider_equilibrated:] == evol[-1]) \
                and step > consider_equilibrated:
                equilibrated = True
                steps_to_equil = np.append(steps_to_equil, [step - consider_equilibrated + 1])                
                print(f"steps to equilibrated: {steps_to_equil[-1]} in simulation {n_sim+1}")
        

            #skip long calculations
            if step > max_steps:
                print(f"simulation {n_sim} did not equilibrate")
                equilibrated = True
            
            if show_anim:
                plt.cla()
                im=plt.imshow(arr, animated=True)
                plt.text(1, -1, f"{method}, step: {step}")
                plt.draw()
                plt.pause(0.01)
            
            
            if method == 'glider':
                indices = np.argwhere(arr == 1)
                com = np.sum(indices, axis = 0)/(len(indices))
                coms = np.append(coms, com)
                            
        if method == 'glider':
            
            coms = coms.reshape(int(len(coms)/2), 2)
            com_mag = np.sqrt(coms[:, 0]**2 + coms[:, 1]**2)
            plt.figure()
            plt.scatter(x = range(len(com_mag)), y = com_mag)
            plt.xlabel("step")
            plt.ylabel("magnitude of com")
            plt.title("Centre of mass")
            plt.show()
                    
            for i in range(len(coms)-1):
                vel = (coms[i+1] - coms[i])
                vels = np.append(vels, vel)
            
            vels = vels.reshape(int(len(vels)/2), 2)
            vels = np.sqrt(vels[:, 0]**2 + vels[:, 1]**2)
            
            
            plt.figure()
            sns.scatterplot(x = range(len(vels)), y = vels)
            plt.xlabel("step")
            plt.ylabel("velocity [cells/step]")
            plt.title("Velocity")
            plt.show()
            
            iso_vel = vels[vels < 1]
            print(f"average velocity: {np.mean(iso_vel)} cells / step")
            
        #write to file
        outfile.write(f"{n_sim} {step}\n")
        
    if show_anim:
        plt.cla()
        im=plt.imshow(arr, animated=True)
        plt.text(1, -1, f"{method}, step: {step}")
        plt.draw()
        plt.show(block=True)
        
    if show_plot:	
        plt.figure()
        sns.histplot(steps_to_equil, bins = 26)
        plt.title("Histogram of steps to equilibrate")
        plt.xlabel("time to steady state")
        plt.ylabel("probability")
        plt.show()

    outfile.close()

main()




end = time.time()
# print(f"time = {end-start}")


#do we need to store data for glider runs

