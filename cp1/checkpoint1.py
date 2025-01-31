import sys
import random
import numpy as np
import pandas as pd
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
    

def full_energy(arr):
    E = 0
    f = arr.reshape(-1)
    for i in range(len(f)):
        for j in range(i+1, len(f)):
            E -= f[i]*f[j]
    return E

def main(auto = False, dynamics_method = "glauber", N = 15, T = 2, arr = None,\
         nsweep = 100, show_anim = True, log_freq = -1):
    
    show_nth = 10
    
    if not auto:
        
        if len(sys.argv) != 4:
            print("%run checkpoint1.py <glauber/kawasaki> <N> <T>")
            return
        
        
        f, dynamics_method, N, T = sys.argv
        N, T = int(N), float(T)
    
    
        arr = set_up(N)

    complete_E = full_energy(arr)
    sweep = N**2
    nsteps = nsweep * sweep
    
    #init fig
    if show_anim:
        fig = plt.figure()
        im=plt.imshow(arr, animated=True)
    
    
    
    if dynamics_method == "glauber":
            
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
            
            complete_E += delta_E
            
            
            if delta_E <= 0:
                arr[site[0], site[1]] *= -1
                
            
            else:
                prob = np.exp(-delta_E / T)
                if random.uniform(0, 1) < prob:
                    arr[site[0], site[1]] *= -1


            #update anim
            if n % show_nth == 0 and show_anim:

                plt.cla()
                im=plt.imshow(arr, animated=True)
                plt.draw()
                plt.pause(0.0001)
            
            #ie store data
            if log_freq > 1 and n/(N*N) % log_freq == 0:
                m1 = np.mean(arr**2)
                m2 = np.mean(arr)**2
                outfile.write(f" {dynamics_method} {T} {m1} {m2} "
                              f"{(m1-m2)/(N*T)} {complete_E/(N**2)} E2 C(E1, E2)\n")                      
                    

    elif dynamics_method == "kawasaki":
       

            
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
                complete_E += net_delta_E
             
            if (np.linalg.norm(np.asarray(sites[0]-sites[1], \
                                          dtype = "float"))) == 1 :
                net_delta_E += 4
            
            
            if net_delta_E <= 0:
                arr[tuple(sites[0])] *= -1
                arr[tuple(sites[1])] *= -1
            
            else:
                prob = np.exp(-net_delta_E / T)
                
                if random.uniform(0, 1) < prob:
                    arr[tuple(sites[0])] *= -1
                    arr[tuple(sites[1])] *= -1


            #update anim
            if n % show_nth == 0 and show_anim:

                plt.cla()
                im=plt.imshow(arr, animated=True)
                plt.draw()
                plt.pause(0.0001)
            
            #ie store data
            if log_freq > 1 and n/(N*N) % log_freq == 0:
                m1 = np.mean(arr**2)
                m2 = np.mean(arr)**2
                outfile.write(f" {dynamics_method} {T} {m1} {m2} {(m1-m2)/(N*T)}"
                              f"{complete_E/(N**2)} E2 C(E1, E2)\n")

    return arr

    


global temps
temps = np.arange(1, 3.1, 0.1)
def collect_data():
    
    N=10
    global outfile
    
    outfile = open('simulation_results.txt', 'w')
    outfile.write('method T <M**2> <M>**2 X <E**2> <E>**2 C\n')
                   
    # 'pre equilibrated'
    arr = np.ones([N, N], dtype=int)
    
    
    
    
    # GLAUBERGLAUBERGLAUBERGLAUBERGLAUBERGLAUBERGLAUBERGLAUBERGLAUBERGLAUBER

    
    #equilibrate
    arr = main(auto = True, dynamics_method = "glauber", N = N, T = 2, \
               arr = arr, nsweep = 100, show_anim = False, log_freq = -1)
    
    
    # 1000 measurements per temperature
    # 10 sweeps in between measurement
    # 1 and 3 in steps of 0.1
    global temps
    temps = np.arange(1, 3.1, 0.1)
    for t in temps:
        "{:.1f}".format(t)
        arr = main(auto = True, dynamics_method = "glauber", N = N, T = t, \
                   arr = arr, nsweep = 100, show_anim = False, log_freq = -1)
        arr = main(auto = True, dynamics_method = "glauber", N = N, T = t, \
                   arr = arr, nsweep = 1000, show_anim = False, log_freq = 10)
                              #10 * num measurements
        print(t)
    outfile.close()
    
collect_data()




def reading_data():
    
    cols = ['method','T','<M**2>','<M>**2', 'X','<E**2>','<E>**2', 'C', 'why']
    data = pd.read_csv('simulation_results.txt', sep=" ", header = 0, names=cols)
    
    df = pd.DataFrame(data)

    print(df[(df['method'] == 'glauber') & (df['T'] == 1)]['X'].mean())


    chi = np.zeros(len(temps))
    for i in range(len(temps)):
        chi[i] = df[(df['method'] == 'glauber') & (df['T'] == temps[i])]['X'].mean()

    plt.scatter(temps, chi)
    plt.show()
    
    
    
    
reading_data()

                  

    




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






#get rid of that line
#work out the visuals


# for kawasaki - consider if probability works - flib both? or check twice?


#change glauber tuple
#change temps loop
#initial equilibration state





