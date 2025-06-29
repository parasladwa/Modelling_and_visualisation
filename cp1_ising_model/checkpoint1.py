import sys
import time
import scipy
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#from astropy.stats import jackknife_stats


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
    
    show_nth = 500

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
    
    FILENAME = '.txt'
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
    
    FILENAMES = ['KAWASAKI_50_1000.txt', 'GLAUBER_50_1000.txt']
    METHODS = ['kawasaki', 'glauber']
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
            
            #
            if METHOD == "glauber":
                print(f" <e2> {E1} <e>2 {E2}")
            #
            C_dict[key] = (E1 - E2) / (N*(key**2))
            
            if METHOD == 'kawasaki':
                C_dict[key] = C_dict[key] / (N**2)
            
        
        
        E_err = {float(t) : float(0) for t in temps}
        for t in temps:
            values = np.array(E_dict[t])
            er1, er2 = np.mean(values**2), np.mean(values)**2
            error = ( (er1-er2) / 1000 )**(1/2)
            E_err[t] = error
        
        
        C_err = jackknife(E_dict, C_dict, 'C')
        
        
        #plt.scatter(C_dict.keys(), C_dict.values())
        plt.errorbar(C_dict.keys(), C_dict.values(), yerr = list(C_err.values()), ecolor = 'r')
        plt.title(f"Specific Heat - {METHOD.capitalize()}")
        plt.show()
        
        
        plt.errorbar(E_average_dict.keys(), E_average_dict.values(),\
                     yerr = list(E_err.values()), ecolor = 'r')
        plt.title(f"Average Energy {METHOD.capitalize()}")
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
            plt.errorbar(X_dict.keys(), X_dict.values(), list(X_err.values()), ecolor = 'r')
            #plt.plot(X_dict.keys(), X_dict.values())
            plt.title(f"Susceptibility - {METHOD.capitalize()}")
            plt.show()
            
            #plt.plot(M_average_dict.keys(), M_average_dict.values())
            plt.errorbar(M_average_dict.keys(), M_average_dict.values(), \
                         yerr = list(M_err.values()), ecolor = 'r')
            plt.title(f"Average Absolute Magnetisation - {METHOD.capitalize()}")
            plt.show()
    
        
        
        
        NET_DATA.append([METHOD, E_average_dict, E_err, C_dict, C_err, M_average_dict, M_err, X_dict, X_err])
        
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
        
plots_and_logging()
    
    
    
    
    
    
    
    
    
    

def create_plots_old():
    
    FILENAME = 'KAWASAKI_50_1000.txt'
    METHOD = 'kawasaki'
    N = 50

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
        
    plt.plot(C_dict.keys(), C_dict.values())
    plt.title(f"Specific Heat - {METHOD.capitalize()}")
    plt.show()
    
    plt.plot(E_average_dict.keys(), E_average_dict.values())
    plt.title(f"Average Energy - {METHOD.capitalize()}")
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
    

    
    
    plt.plot(X_dict.keys(), X_dict.values())
    plt.title(f"Susceptibility - {METHOD.capitalize()}")
    plt.show()
    
    plt.plot(M_average_dict.keys(), M_average_dict.values())
    plt.title("Average Absolute Magnetisation - {METHOD.capitalize()}")
    plt.show()

#create_plots()    





end = time.time()
print(f"net time ; {end-start}")
                  





def reading_data_working_old():
    
    #read file
    with open('GLAUBER_50_1000.txt', 'r') as f:
        lines = f.readlines()
    data = []
    for line in lines:
        data.append(line.split())
        
    
    #plot of the susceptibility and of the average absolute value of the magnetisation
    M_dict_g = {float(t): [] for t in temps}

    for line in data:
        if line[0] == "glauber":
            M_dict_g[float(line[1])].append(float(line[2]))
        
    
    X_dict_g = {float(t): [] for t in temps}
    M_average_dict_g = {float(t): [] for t in temps}
    for key in M_dict_g.keys():
        
        Ms_g = np.array(M_dict_g[key])    
        M_average_dict_g[key] = abs(np.mean(Ms_g))
        
        m1, m2 = np.mean(Ms_g**2), np.mean(Ms_g)**2
        N = 50
        X_dict_g[key] = (m1-m2) / (N*key)
    
    plt.plot(X_dict_g.keys(), X_dict_g.values())
    plt.title('Glauber - Susceptibility')
    plt.show()

    plt.plot(M_average_dict_g.keys(), M_average_dict_g.values())
    plt.title('Glauber - Average Absolute Value of the Magnetisation')
    plt.show()
    
    
    #plot of specific heat and average energy
    E_dict_g = {float(t): [] for t in temps}
    for line in data:
        if line[0] == "glauber":
            E_dict_g[float(line[1])].append(float(line[3]))
    
    C_dict_g = {float(t): [] for t in temps}
    E_average_dict_g = {float(t): [] for t in temps}
    for key in E_dict_g:
        
        Es_g = np.array(E_dict_g[key])
        E_average_dict_g[key] = np.mean(Es_g)
        
        E1, E2 = np.mean(Es_g**2), np.mean(Es_g)**2
        C_dict_g[key] = (E1 - E2) / (N*(key**2))
   

    plt.plot(C_dict_g.keys(), C_dict_g.values())
    plt.title("Specific Heat")
    plt.show()
    
    plt.plot(E_average_dict_g.keys(), E_average_dict_g.values())
    plt.title("Glauber - Average Energy")
    plt.show()
    

def reading_data1():
    
    cols = ['method','T','<M**2>','<M>**2', 'X','<E**2>','<E>**2', 'C', 'why']
    data = pd.read_csv('simulation_results2.txt', sep=" ", header = 0, names=cols)
    
    df = pd.DataFrame(data)

    print(df[(df['method'] == 'glauber') & (df['T'] == 1)]['X'].mean())


    chi = np.zeros(len(temps))
    for i in range(len(temps)):
        chi[i] = df[(df['method'] == 'glauber') & (df['T'] == temps[i])]['X'].mean()

    plt.scatter(temps, chi)
    plt.show()


def reading_data2():
    
    with open('total.txt', 'r') as f:
        lines = f.readlines()
    
    data = []
    for l in lines:
        p = l.split()
        data.append(p)
        
        
    print(data[1])
    
    
    
    energy_values = {float(t): [] for t in temps}
    for line in data:
        if line[0] == "glauber":
            energy_values[float(line[1])].append(float(line[5]))
    
    
    N=50
    C_values = {float(t): 0 for t in temps}
    for t in temps:
        E_arr = np.asarray(energy_values[t], dtype = 'float')
        E1 = np.mean(E_arr**2)
        E2 = np.mean(E_arr)**2
        C = (E1 - E2) / (N* (t**2))
        
        C_values[t] = (C)
    
    plt.plot(C_values.keys(), C_values.values())
    plt.title("specific heat capacity")
    plt.show()


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




# THAT ONE GRAPH
# PRINT ERRORS
#
#
#
#
#
#
#
#
