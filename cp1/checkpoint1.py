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
        print("%run checkpoint1.py <glauber/kawasaki> <N> <T>")
        return
    
    
    f, dynamics_method, N, T = sys.argv
    N, T = int(N), float(T)
    
    arr = set_up(N)
    #init fig
    fig = plt.figure()
    im=plt.imshow(arr, animated=True)
    
    
    
    if dynamics_method == "glauber":
        
        #list for multiple
        #for glauber want to use 1<T<3, dt = 0.1
        #use just one for now
        temps = np.array([T])
        
        for t in temps:
            
            nsteps=250000######################################HOW MANY
            
            for n in range(nsteps):
                
                #choose site randomly
                site = np.array([random.randint(0, N-1), \
                                 random.randint(0, N-1)])
                
                #find neighbours addresses
                nn_addresses = find_nn(site, N)
    
                #now find spins at nn's
                neighbours = np.zeros([4])
                for i in range(len(nn_addresses)):
                    neighbours[i] = arr[int(nn_addresses[i][0])] \
                                        [int(nn_addresses[i][1])]
                
                #find energy & energy change
                #change = after - initial
                energy_init = find_energy(neighbours, arr, site)
                energy_fin = -1 * energy_init
                delta_E = energy_fin - energy_init
                
                
                """
                branch here
                
                case1:
                    delta E <= 0
                    flip occurs
                    
                case2:
                    flip with probability of exp(-delta E / k_b T)
                """
                
                
                if delta_E <= 0:
                    arr[site[0], site[1]] *= -1
                    
                
                else:
                    prob = np.exp(-delta_E / t)
                    if random.uniform(0, 1) < prob:
                        arr[site[0], site[1]] *= -1


                #update anim
                if n%10 == 0:

                    plt.cla()
                    im=plt.imshow(arr, animated=True)
                    plt.draw()
                    plt.pause(0.0001)
                
                    
                    M = np.sum(arr)
                    print(M)
                    



    elif dynamics_method == "kawasaki":
       
        
        temps = np.array([T])
        
        for t in temps:
            
            nsteps = 2500000
            
            for n in range(nsteps):
                
                #choose 2 sites randomly
                sites = np.array([[random.randint(0, N-1), \
                                   random.randint(0, N-1)], \
                                  [random.randint(0, N-1), \
                                   random.randint(0, N-1)]])
                    
                #eliminate equal spins case
                if arr[tuple(sites[0])] == arr[tuple(sites[1])]:
                    continue
                
                
                #consider consecutive spinflips
                net_delta_E = 0
                for site in sites:
                    
                    neighbours = np.zeros([4])
                    nn_addresses = find_nn(site, N)
                    
                    for i in range(len(nn_addresses)):
                        neighbours[i] = arr[int(nn_addresses[i][0])] \
                                            [int(nn_addresses[i][1])]
                            

                    #find energy & energy change
                    #change = after - initial
                    energy_init = find_energy(neighbours, arr, site)
                    energy_fin = -1 * energy_init
                    delta_E = energy_fin - energy_init
                    
                    net_delta_E += delta_E
                
                
                if net_delta_E <= 0:
                    arr[tuple(sites[0])] *= -1
                    arr[tuple(sites[1])] *= -1
                
                
                else:
                    prob = np.exp(-delta_E / t)
                    
                    if random.uniform(0, 1) < prob:
                        arr[tuple(sites[0])] *= -1
                        arr[tuple(sites[1])] *= -1



                #update anim
                if n%500 == 0:

                    plt.cla()
                    im=plt.imshow(arr, animated=True)
                    plt.draw()
                    plt.pause(0.0001)
                    
            
                    
                    M = np.sum(arr)
                    print(M)
            
            
            
            
                    
def testing_funcs(N):
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
    energy_init = find_energy(neighbours, arr, site)
    
    # change = after - initial
    energy_fin = -1 * energy_init
    delta_E = energy_fin - energy_init
    print(delta_E)
    
    print(arr)
    
    
    """
    branch here
    
    case1:
        delta E <= 0
        flip occurs
        
    case2:
        flip with probability of exp(-delta E / k_b T)
    """
    
    if delta_E <= 0:
        arr[site[0], site[1]] *= -1
        print("E loss flipped")
    
    else:
        prob = np.exp(-delta_E / T)
        print(prob)
        
        #consider flip (rng)
        if random.uniform(0, 1) < prob:
            arr[site[0], site[1]] *= -1
            print("rng flipped")


def more_testing(N):
    

    arr = set_up(N)

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
    
    
    print(delta_E)
    if delta_E > 0:
        t= 2
        prob = np.exp(-delta_E / t)
        print(prob)
        p=random.uniform(0, 1)
        if p < prob:
            print(p, prob)
            print('flipppppeeeeeeedddddddddddddddddddd')
            arr[site[0], site[1]] *= -1


if __name__ == "__main__":
    main()




#get rid of that line
#work out the visuals
#sweeps????????????????????


# for kawasaki - consider if probability works - flib both? or check twice?
# must it be command line ?


# same function or too much hastle?
# handle false inputs 'glaber' ?

#change glauber tuple





