import time
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap











def initialise(N):
    return np.random.choice([0, 1, 2], size=(N, N))












def simulate(arr, p1, p2, p3, N, nsweeps, show_anim, show_nth, case, log, fractions = False):

    if fractions:
        fractions_filename = 'data_fractions.txt'
        outfile_fractions = open(fractions_filename, 'w')
        outfile_fractions.write("<I> <S> <R>\n")
    
    
    nsteps = nsweeps*N*N

    average_infected = np.zeros(nsweeps)
    average_susceptible = np.zeros(nsweeps)
    average_recovered = np.zeros(nsweeps)
    
    if show_anim:
        cmap = ListedColormap(["white", "darkred", "grey"])
        fig = plt.figure()
        im=plt.imshow(arr, animated=True)

        
    for step in range(nsteps):

        i, j = np.random.randint(0, N, 2)
        
        # S0 -> I1
        if arr[i, j] == 0:
            neighbors = [
                ((i - 1) % N, j),
                ((i + 1) % N, j),
                (i, (j - 1) % N),
                (i, (j + 1) % N) 
            ]

            if any(arr[x, y] == 1 for x, y in neighbors):
                if np.random.rand() < p1:
                    arr[i, j] = 1

        # I1 -> R2
        elif arr[i, j] == 1:
            if np.random.rand() < p2:
                arr[i, j] = 2
        
        # R2 -> S0
        else:
            if np.random.rand() < p3:
                arr[i, j] = 0
    

        if log and step % (N*N) == 0:
            average_infected[step//(N*N)] = np.sum(arr == 1)
            average_susceptible[step//(N*N)] = np.sum(arr == 0)
            average_recovered[step//(N*N)] = np.sum(arr == 2)
        
        
        if step % show_nth == 0:
            
            if fractions:
                outfile_fractions.write(f"{np.sum(arr == 1)} {np.sum(arr == 0)} {np.sum(arr == 2)}\n")
            
            if show_anim:
                plt.cla()
                im=plt.imshow(arr, animated=True, cmap=cmap)
                plt.title(f"SIRS Model, {case} case\nsweep: {step/(N*N)}, p1: {p1}, p2: {p2}, p3: {p3}")
                plt.draw()
                plt.pause(0.01)
            

    
    
    return arr, average_infected, average_susceptible, average_recovered














def phase_plane_plot():
    
    filename = "test.txt"
    
    with open(filename, 'r') as file:
        lines = file.readlines()
       
    data = []

    for l in lines:
        l = l.strip()
        l = l.split()
        data.append(l)
        
    data = np.array(data)

    p1s = data[1:, 0]
    p2s = data[1:, 1]
    p3s = data[1:, 2]
    Is = data[1:, 3]
    I2s = data[1:, 4]
    

    all_ps = np.arange(0, 1.05, 0.05, dtype=float)

    
    mapped = np.zeros((len(all_ps), len(all_ps)), dtype= float)


    for dp in data[1:]:
        
        dp = np.array(dp, dtype= float)
        
        p1, p3 = dp[0], dp[2]

        
        i = np.where( all_ps == p1)
        j = np.where( all_ps == p3)
        
        mapped[i[0][0], j[0][0]] = dp[3] / 2500
    
    
    labels = np.round(all_ps, decimals=2)
    
    plt.figure()
    plt.title('<I>/N for p1 p3 plane')
    sns.heatmap(mapped, xticklabels=labels, yticklabels=labels)
    plt.gca().invert_yaxis()
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
def variance_plane_plot():
    
    filename = "test.txt"
    
    with open(filename, 'r') as file:
        lines = file.readlines()
       
    data = []

    for l in lines:
        l = l.strip()
        l = l.split()
        data.append(l)
        
    data = np.array(data)

    p1s = data[1:, 0]
    p2s = data[1:, 1]
    p3s = data[1:, 2]
    Is = data[1:, 3]
    I2s = data[1:, 4]
    

    all_ps = np.arange(0, 1.05, 0.05, dtype=float)

    
    mapped = np.zeros((len(all_ps), len(all_ps)), dtype= float)


    for dp in data[1:]:
        
        dp = np.array(dp, dtype= float)
        
        p1, p3 = dp[0], dp[2]

        
        i = np.where( all_ps == p1)
        j = np.where( all_ps == p3)
        
        mapped[i[0][0], j[0][0]] = (1/2500) * (dp[4] - dp[3]**2)
    
    
    labels = np.round(all_ps, decimals=2)
    
    plt.figure()
    sns.heatmap(mapped, xticklabels=labels, yticklabels=labels)
    plt.title('variance in I for p1 p3 plane')
    plt.gca().invert_yaxis()
    plt.show()
    
    
    
    
    
    
    
    
    

def plot_cut():
    
    f_name = 'data_from_cut.txt'
    data = np.loadtxt(f_name, comments = 'p', dtype=float)

    p1s = data[:, 0]
    p2s = data[:, 1]
    p3s = data[:, 2]

    I = data[:, 3]
    Isqr = data[:, 4]
    
    vals = np.zeros(len(I))
    
    for i in range(len(vals)):
        vals[i] = (1/2500)*(Isqr[i] - I[i]**2)
    
    plt.figure()
    plt.plot(p1s, vals)
    plt.title('variance plot for fixed p2 p3')
    plt.xlabel('p1')
    plt.ylabel('variance')
    plt.show()





def fractions_plots():
   
    f_name = 'data_fractions.txt'
    data = np.loadtxt(f_name, comments = '<', dtype=float)
    data /= 2500
    Is = data[1:, 0]
    Ss = data[1:, 1]
    Rs = data[1:, 2]

    x = np.array(list(range(len(Is))))
    
    plt.figure()
    plt.scatter(x, Is, c='r', s=1.5)
    plt.scatter(x, Ss, c='black', s=1.5)
    plt.scatter(x, Rs, c='grey', s=1.5)
    plt.xlabel('steps but different')
    plt.ylabel('fraction')
    
    plt.title("fraction of states")
    plt.text(x[-1]*1.07, 0.55,"S", c = 'black')
    plt.text(x[-1]*1.07, 0.5,"I", c = 'r')
    plt.text(x[-1]*1.07, 0.45,"R", c = 'grey')
    plt.show()
    
    
    
    







def main():
    # S0 I1 R2
    
    cases = {
        'absorbing' : np.array([0.1, 0.9, 0.9], dtype=np.float64),
        'dynamic_eq' : np.array([0.7, 0.5, 0.5], dtype=np.float64),
        'cyclic' : np.array([0.8, 0.1, 0.01], dtype=np.float64)
    }
    
    parser = argparse.ArgumentParser(description="SIR Model Simulation with CLI Inputs")
    parser.add_argument("-a", "--auto", action='store_true', help="full auto mode for data collection")
    parser.add_argument("-ac", "--auto_cut", action = 'store_true', help = "auto mode for data collection along specified cut")
    parser.add_argument("-N", "--N", type=int, default=50, help="Size of the lattice")
    parser.add_argument("-s", "--show_anim", action= "store_true", help="Show animation of the simulation")
    parser.add_argument("-ns", "--nsweeps", type=int, default=5, help="Number of sweeps")
    parser.add_argument("-p1", "--p1", type=float, default=0.5, help="Probability of S -> I transition")
    parser.add_argument("-p2", "--p2", type=float, default = 0.5, help="Probability of I -> R transition")
    parser.add_argument("-p3", "--p3", type=float, default = 0.5, help="Probability of R -> S transition")
    parser.add_argument("-c", "--case", choices=["absorbing", "dynamic_eq", "cyclic"], help="Choose a predefined case")
    parser.add_argument("-nth", "--show_nth", type=int, default=2500, help="Show every nth step")
    parser.add_argument("-phase", "--p1_p3_phase", action = 'store_true', help = "Plots <I> / N phase across p1, p3 plane")
    parser.add_argument("-var", "--phase_variance", action = 'store_true', help = "Plots ( <I2> - <I>2 ) / N phase across p1, p3 plane")
    parser.add_argument("-f", "--fractions", action = 'store_true', help = "Plots the fractions of each state across the simulation")
    parser.add_argument("-pc", "--plot_cut", action = 'store_true', help = "Plots the variance of cut simulation")



    args = parser.parse_args()

    if args.plot_cut:
        plot_cut()
        return
    
    
    if args.p1_p3_phase:
        phase_plane_plot()
        return
    
    
    if args.phase_variance:
        variance_plane_plot()
        return
    
    if args.auto_cut:
        
        p1s = np.arange(0.2, 0.51, 0.01)
        p2s = np.array([0.5], dtype=float)
        p3s = np.array([0.5], dtype=float)

        N = 50
        show_anim = args.show_anim
        equilibration_sweeps = 100
        measurement_sweeps = 10
        
        #unused
        show_nth = 100
        case = 'Default'
        
        cut_filename = 'data_from_cut.txt'
        outfile_cut = open(cut_filename, 'w')
        outfile_cut.write("p1 p2 p3 <I> <I^2> <S> <R>\n")

    elif args.auto:
        
        # all possible ps
        all_ps = np.arange(0, 1.05, 0.05)
        p2s = np.array([0.5], dtype=float)
        p1s = all_ps
        p3s = all_ps
        
        N = 50
        show_anim = args.show_anim
        equilibration_sweeps = 10
        measurement_sweeps = 10
        
        #unused
        show_nth = 100
        case = 'Default'
  
        filename = "test.txt"
        outfile_plane = open(filename, 'w')
        outfile_plane.write("p1 p2 p3 <I> <I^2> <S> <R>\n")
    
    
    
    

    else:
        if args.case == None:
            case = 'Default'
            p1, p2, p3 = args.p1, args.p2, args.p3
        else: 
            case = args.case
            p1, p2, p3 = cases[args.case]
        
        N = args.N
        show_anim = args.show_anim
        nsweeps = args.nsweeps
        nsteps = nsweeps*N*N
        show_nth = args.show_nth

        p1s = np.array([p1], dtype = float)
        p2s = np.array([p2], dtype = float)
        p3s = np.array([p3], dtype = float)
    

    for i, p1 in enumerate(p1s):

        for j, p3 in enumerate(p3s):
            
            for k, p2 in enumerate(p2s):
            
            
                arr = initialise(N)
                log = False
                start = time.time()
                
                if args.auto or args.auto_cut:
                    nsweeps = equilibration_sweeps
                    arr, ignore, ignore, ignore= simulate(arr, p1, p2, p3, N, nsweeps, show_anim, show_nth, case, log=False)
                    log = True
                    nsweeps = measurement_sweeps
                
                
                if not args.auto:
                    print(f"simulation started with {case} case "
                    f"p1: {p1}, p2: {p2}, p3: {p3}\n"
                    f"Number of sweeps: {nsweeps}\n\n")
                    
                arr, av_infected, average_susceptible, average_recovered = simulate(arr, p1, p2, p3, N, nsweeps, show_anim, show_nth, case, log, fractions = args.fractions)
                end = time.time()
                print(f"time = {end-start}")
                
                if args.auto:
                    outfile_plane.write(f"{p1} 0.5 {p3} {np.mean(av_infected)} {np.mean(av_infected**2)} {np.mean(average_susceptible)} {np.mean(average_recovered)}\n")
                    print(f"{p1} 0.5 {p3} {np.mean(av_infected)} {np.mean(av_infected**2)} {np.mean(average_susceptible)} {np.mean(average_recovered)}\n")
                
                if args.auto_cut:
                    outfile_cut.write(f"{p1} 0.5 {p3} {np.mean(av_infected)} {np.mean(av_infected**2)} {np.mean(average_susceptible)} {np.mean(average_recovered)}\n")
                    print(f"{p1} 0.5 {p3} {np.mean(av_infected)} {np.mean(av_infected**2)} {np.mean(average_susceptible)} {np.mean(average_recovered)}\n")


    
    if args.fractions:
        fractions_plots()
main()


#close outfiles
#label axes
