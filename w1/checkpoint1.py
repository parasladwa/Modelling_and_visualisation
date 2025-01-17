import sys
import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

N = 50


def set_up():
    #randomly set up array
    arr = np.zeros([N, N], dtype=int)
    
    for i in range(0, len(arr)):
        for j in range(0, len(arr[0])):
            arr[i][j] = random.choice([-1, 1])
    
    return arr
    


def find_nn(site):
    
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
    '''
    so far this is to test functions
    will be in a loop later once funcs are completed

    '''
    #random arr
    arr = set_up()    
    
    #choose site randomly
    site = [random.randint(0, N-1), random.randint(0, N-1)]
    site = np.array([0, 0], dtype = int)
    
    #find neighbours POOR VARIABLE NAME , THIS IS NEIGHBOUR ADRESSES NOT NEIGHBOURS fix this
    nn_addresses = find_nn(site)
    
    
    #now find values at nn's
    neighbours = np.zeros([4])
    for i in range(len(nn_addresses)):
        neighbours[i] = arr[int(nn_addresses[i][0])][int(nn_addresses[i][1])]
    
    
    #find energy & energy change
    energy = find_energy(neighbours, arr, site)
    print(energy, 'energy')
        
    

    
    
    #find energy of site
    energy = 0
    
    
    





main()