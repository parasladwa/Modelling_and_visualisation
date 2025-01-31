import sys
import random
import numpy as np
import matplotlib.pyplot as plt

def set_up(N):
    #randomly set up array
    arr = np.zeros([N, N], dtype=int)
    
    for i in range(0, len(arr)):
        for j in range(0, len(arr[0])):
            arr[i][j] = random.choice([-1, 1])
    
    return arr
    


def find_nn(site, N):
    
    #find addresses of nearest neighbors
    neighbours = np.zeros(shape=(4,2))
    
    #down
    neighbours[0] = site + np.array([1, 0])
    #up
    neighbours[1] = site + np.array([-1, 0])
    #right
    neighbours[2] = site + np.array([0, 1])
    #left
    neighbours[3] = site + np.array([0, -1])
    
    #boundary conditions
    for i in range(len(neighbours)):
        for j in range(len(neighbours[0])):
            
            if neighbours[i][j] == N:
                neighbours[i][j] = 0
                
            if neighbours[i][j] == -1:
                neighbours[i][j] = N-1
                
            neighbours[i][j] = np.asarray(neighbours[i][j], dtype=int) # check if i can delete this
            
    return neighbours


    
def find_energy(nn, arr, site):
    
    
    
    
    #find energy given nearest neighbours (np.arr)
    
    # E(S_site) = - J * sum [ S_site * S_nn ]
    # J = 1
    site_spin = arr[site[0], site[1]]
    energy = 0
    
    for spin in nn:
        energy -= site_spin * spin
        
    return energy




def main():
    
    
    if len(sys.argv) != 4:
        print("%run checkpoint1part2.py <glauber/kawasaki> <N> <T or auto>")
        return
    
    f, dynamics_method, N, T = sys.argv
    N = int(N)
    
    
    if T == 'auto':
        temps = np.arange(1, 3.1, 0.1)
        temps = np.arange(1, 1.3, 0.1)
        #why are they slightly off
        
    else:
        temps = np.array([T])
    
    
    for t in temps:
        t = round(t, 2)
        
        
        nsteps = 2500
        for n in range(nsteps):
            
            #choose site randomly
            site = [random.randint(0, N-1), random.randint(0, N-1)]
    
            #find neighbours addresses
            nn_addresses = find_nn(site, N)
            
            #now find spins at nn's
            neighbours = np.zeros([4])
            for i in range(len(nn_addresses)):
                neighbours[i] = arr[int(nn_addresses[i][0])][int(nn_addresses[i][1])]
            
            #find energy & energy change
            #change = after - initial
            energy_init = find_energy(neighbours, arr, site)
            energy_fin = -1 * energy_init
            delta_E = energy_fin - energy_init
            
            
            
            if dynamics_method == "kawasaki":
                #compute a flip again
                
            
            
    
main()