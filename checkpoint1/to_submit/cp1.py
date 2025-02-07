import sys
import time
import scipy
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.stats import jackknife_stats


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
         nsweep = 1000, show_anim = True, log_freq = -1):
    
    show_nth = 10*2500

    if not auto:
        
        if len(sys.argv) != 4:
            print("%run cp1.py <glauber/kawasaki> <N> <T>")
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
            
            
            
            
            if delta_E <= 0:
                arr[site[0], site[1]] *= -1
                complete_E += delta_E
            
            else:
                prob = np.exp(-delta_E / T)
                if random.uniform(0, 1) < prob:
                    arr[site[0], site[1]] *= -1
                    complete_E += delta_E


            #update anim
            if n % show_nth == 0 and show_anim:

                plt.cla()
                im=plt.imshow(arr, animated=True)
                plt.draw()
                plt.pause(0.0001)
            
            #ie store data
            if log_freq > 1 and n/(N*N) % log_freq == 0:
                
                M = np.sum(arr)
                outfile.write(f" {dynamics_method} {T} {M} "
                              f"{complete_E}\n")                      
                   
                
                
                
                
    
    if dynamics_method == "kawasaki":
                   

        for n in range(nsteps):

            
            #choose 2 sites randomly
            sites = np.array([[random.randint(0, N-1), \
                               random.randint(0, N-1)], \
                              [random.randint(0, N-1), \
                               random.randint(0, N-1)]])
                
            #eliminate equal spins case
            if arr[tuple(sites[0])] == arr[tuple(sites[1])]:
                continue
           
            
            
            # E before
            nn_addresses0 = find_nn(sites[0], N)
            neighbours0 = np.zeros([4])
            for i in range(len(nn_addresses0)):
                neighbours0[i] = arr[int(nn_addresses0[i][0])] \
                                    [int(nn_addresses0[i][1])]
            nn_addresses1 = find_nn(sites[1], N)
            neighbours1 = np.zeros([4])
            for i in range(len(nn_addresses1)):
                neighbours1[i] = arr[int(nn_addresses1[i][0])] \
                                    [int(nn_addresses1[i][1])]
            energy_init = find_energy(neighbours0, arr, sites[0]) + \
                            find_energy(neighbours1, arr, sites[1])
             
                
            #swap and E after
            arr[tuple(sites[0])], arr[tuple(sites[1])] = arr[tuple(sites[1])], arr[tuple(sites[0])]

            
            neighbours0 = np.zeros([4])
            for i in range(len(nn_addresses0)):
                neighbours0[i] = arr[int(nn_addresses0[i][0])] \
                                    [int(nn_addresses0[i][1])]
            neighbours1 = np.zeros([4])
            for i in range(len(nn_addresses1)):
                neighbours1[i] = arr[int(nn_addresses1[i][0])] \
                                    [int(nn_addresses1[i][1])]

            energy_final = find_energy(neighbours0, arr, sites[0]) + \
                            find_energy(neighbours1, arr, sites[1])
            
            
            delta_E = energy_final - energy_init
            
            
            
            arr[tuple(sites[0])], arr[tuple(sites[1])] = arr[tuple(sites[1])], arr[tuple(sites[0])]

            if delta_E <= 0:
                complete_E += delta_E
                arr[tuple(sites[0])], arr[tuple(sites[1])] = arr[tuple(sites[1])], arr[tuple(sites[0])]
            
            else:
                prob = np.exp(-delta_E / T)
               
                if random.uniform(0, 1) < prob:
                    complete_E += delta_E
                    arr[tuple(sites[0])], arr[tuple(sites[1])] = arr[tuple(sites[1])], arr[tuple(sites[0])]
            
                    


            #update anim
            if n % show_nth == 0 and show_anim:
                plt.cla()
                im=plt.imshow(arr, animated=True)
                plt.draw()
                plt.pause(0.0001)
            
            #ie store data
            if log_freq > 1 and n/(N*N) % log_freq == 0:
                
                M = np.sum(arr)
                outfile.write(f" {dynamics_method} {T} {M} "
                              f"{complete_E}\n")
    return arr

main()





start = time.time()
global temps
temps = np.arange(1, 3.1, 0.1)




def collect_data():

    
    global outfile
    
    FILENAME = 'kawasaki_50_1000.txt'
    METHOD = 'kawasaki'
    MEASUREMENTS = 1000
    N=50
    
    

    confirm = input(f"      method : {METHOD}\n"
                    f"           N : {N}\n"
                    f"measurements : {MEASUREMENTS}\n"
                    f"    filename : {FILENAME}\n\n"
                    f"CONFIRM ??? [y/n] :")
    
    if confirm != 'y':
        print('terminated')
        return
    
    
    if METHOD == 'glauber':
        arr = np.ones((N, N))
    elif METHOD == 'kawasaki':
        arr = set_up(N)
    else:
        print('terminated - incorrect method')
        return
    
    
    outfile = open(FILENAME, 'w')
    outfile.write('method T M E\n')

    # 10 sweeps in between measurement
    # 1 and 3 in steps of 0.1
    global temps
    temps = np.arange(1, 3.1, 0.1)
    for t in temps:
        print(t)
        start1= time.time()
        arr = main(auto = True, dynamics_method = METHOD, N = N, T = t, \
                   arr = arr, nsweep = 100, show_anim = False, log_freq = -1)
        arr = main(auto = True, dynamics_method = METHOD, N = N, T = t, \
                   arr = arr, nsweep = 10 * MEASUREMENTS, show_anim = False,\
                       log_freq = 10)
        
        end1 = time.time()
        print(f"time ; {end1-start1}")
    print(t)
    
    outfile.close()
    
    
#collect_data()










def jackknife(data, true_data, which):
    C_err = {float(t):float(0) for t in data.keys()}
    for t in data.keys():
        
        
        temp_c = data[t]
        Cs = np.array([], dtype = float)
        true = true_data[t]
        
        for i in range(len(temp_c)):
            temp_data = np.delete(temp_c, i)
            
            
            e1, e2 = np.mean(temp_data**2), np.mean(temp_data)**2
            if which == 'C':
                C = (e1 - e2) / (50 * (t**2))
            if which == 'X':
                C = (e1 - e2) / (50 * t)
            Cs = np.append(Cs, C)
                
        s = 0
        for ci in Cs:
            s += (ci - true)**2
        C_err[t] = (s)**(1/2)
        
        
    return C_err
        
        
            




def plots_and_logging():
    
    FILENAMES = ['glauber_50_1000.txt', 'kawasaki_50_1000.txt']
    METHODS = ['glauber', 'kawasaki']
    N = 50
    
    
    outfile_name = "all_results.txt"
    datafile = open(outfile_name, 'w')
    datafile.write('method T E E_err C C_err M M_err X X_err\n')

    NET_DATA = []
    for i in range(len(FILENAMES)):
        
        FILENAME, METHOD = FILENAMES[i], METHODS[i]
    
        with open(FILENAME, 'r') as f:
            lines = f.readlines()
        
        data = []
        for line in lines:
            data.append(line.split())
        
            
    
        
        #specific heat and average energy
        E_dict = {float(t):[] for t in temps}
        for line in data:

            if line[0] == METHOD:
                E_dict[float(line[1])].append(float(line[3]))

    
        C_dict = {float(t) : [] for t in temps}
        E_average_dict = {float(t) : None for t in temps} 
        
        for key in E_dict:             
            Es = np.array(E_dict[key])
            E_average_dict[key] = np.mean(Es)
            
            E1, E2 = np.mean(Es**2), np.mean(Es)**2
            
            C_dict[key] = (E1 - E2) / (N*(key**2))
            
            # if METHOD == 'kawasaki':
            #     C_dict[key] = C_dict[key] / (N**2)
            
        
        
        E_err = {float(t) : float(0) for t in temps}
        for t in temps:
            values = np.array(E_dict[t])
            er1, er2 = np.mean(values**2), np.mean(values)**2
            error = ( (er1-er2) / 1000 )**(1/2)
            E_err[t] = error
        
        
        C_err = jackknife(E_dict, C_dict, 'C')
        
        

        plt.errorbar(C_dict.keys(), C_dict.values(),\
                     yerr = list(C_err.values()), ecolor = 'r')
        plt.title(f"Specific Heat - {METHOD.capitalize()}")
        plt.xlabel("Temperature")
        plt.ylabel("Specific Heat Capacity")
        plt.show()
        
        
        plt.errorbar(E_average_dict.keys(), E_average_dict.values(),\
                     yerr = list(E_err.values()), ecolor = 'r')
        plt.title(f"Average Energy {METHOD.capitalize()}")
        plt.xlabel("Temperature")
        plt.ylabel("Average Energy")
        plt.show()
        
        
        
        
        
        
        #susceptibility and average absolute magnetisation

        M_dict = {float(t) : [] for t in temps}
        for line in data:
            if line[0] == METHOD:
                M_dict[float(line[1])].append(float(line[2]))
        
        
        X_dict = {float(t) : [] for t in temps}
        M_average_dict = {float(t) : None for t in temps}
        for key in M_dict:
            
            Ms = np.array(M_dict[key])
            M_average_dict[key] = abs(np.mean(Ms))
            
            M1, M2 = np.mean(Ms**2), np.mean(Ms)**2
            X_dict[key] = (M1 - M2) / (N * key)
        
        
        # M_err
        M_err = {float(t) : float(0) for t in temps}
        for t in temps:
            values = np.array( M_dict[t])
            er1, er2 = np.mean(values**2), np.mean(values)**2
            error = ( (er1-er2) / 1000 )**(1/2)
            M_err[t] = error
        
        
        
        X_err = jackknife(M_dict, X_dict, 'X')
        if METHOD == 'glauber':
            plt.errorbar(X_dict.keys(), X_dict.values(), list(X_err.values()),\
                         ecolor = 'r')
            plt.title(f"Susceptibility - {METHOD.capitalize()}")
            plt.xlabel("Temperature")
            plt.ylabel("Susceptibility")
            plt.show()
            
            plt.errorbar(M_average_dict.keys(), M_average_dict.values(), \
                         yerr = list(M_err.values()), ecolor = 'r')
            plt.title(f"Average Absolute Magnetisation - {METHOD.capitalize()}")
            plt.xlabel("Temperature")
            plt.ylabel("Absolute Magnetisation")
            plt.show()
    
        
        
        
        NET_DATA.append([METHOD, E_average_dict, E_err, C_dict, C_err,\
                         M_average_dict, M_err, X_dict, X_err])
        
        f.close()

    

    for i, d in enumerate(NET_DATA):
        for t in d[1].keys():
            #('method T E E_err C C_err M M_err X X_err\n')
            datafile.write(f"{d[0]} {t} "
                           f"{d[1][t]} {d[2][t]} " # E delE
                           f"{d[3][t]} {d[4][t]} " # C delC
                           f"{d[5][t]} {d[6][t]} " # M delM
                           f"{d[7][t]} {d[8][t]}\n") #X X_err
    
 
        
    datafile.close()
        
#plots_and_logging()
    
    
    
    
    
