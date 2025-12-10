import time
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


# FOR DEMONSTRATION

# -s -ns 300 -f -p1 0.1 -p2 0.2 -p3 0.3

# -s -ns 600 -f -c cyclic
# immunity -pi
# -phase
# -var
# -pc


# 4, 8, 13 is a good seed for cyclic case
SEED = 13
np.random.seed(SEED)






def initialise(N, f_immune):
    """
    args:
        N: int ; size of lattice
        f_immune: float ; fraction of immune sites

    returns:
        arr : np.ndarray
        N by N grid of sites with f_immune*N*N sites
        and randomly initilaised sites of state S I R 
        or 0, 1, 2 respectively

    """
    arr =  np.random.choice([0, 1, 2], size=(N, N))
    if f_immune == 0:
        return arr

    
    n_sites = round(N*N*f_immune)
    
    random_sites = []
    while len(random_sites) < n_sites:
        
        site =  list(np.random.randint(0, N, 2))
        if site not in random_sites:
            random_sites.append(site)
    
    
    for s in random_sites:
        arr[tuple(s)] = 4
    
    return arr








def jackknife(data):
    """
    args:
        data: np.ndarray ; input array of data

    returns:
        float ; jackknife estimate of the standard error 
        in the variance-like quantity (1/2500) * (⟨x²⟩ - ⟨x⟩²)
    
    description:
        Computes a jackknife error estimate for the sample statistic 
        (1/2500) * (mean(data**2) - mean(data)**2). 
        The function systematically removes one data point at a time 
        to generate pseudo-samples, evaluates the statistic for each, 
        and returns the square root of the summed squared deviations 
        from the full-sample value.
    """

    
    true = (1/(2500) * (np.mean(data**2) - (np.mean(data)**2)))
    vals = []
    
    for i in range(len(data)):
        d = np.delete(data, i)
        vals.append((1/2500) * (np.mean(d**2) - (np.mean(d)**2)))
    
    vals = np.array(vals, dtype=float)
    err = 0
    for v in vals:
        err += (v-true)**(2)
    return (err)**(1/2)
        
        
        
        
        









def simulate(arr, p1, p2, p3, N, nsweeps, show_anim, show_nth, case, log, fractions = False):
    """
    args:
        arr: np.ndarray ; initial N x N lattice of site states 
            (0: susceptible, 1: infected, 2: recovered, 4: immune)
        p1: float ; infection probability (S → I transition)
        p2: float ; recovery probability (I → R transition)
        p3: float ; loss-of-immunity probability (R → S transition)
        N: int ; lattice size
        nsweeps: int ; number of full-lattice sweeps to simulate 
                      (1 sweep = N*N update attempts)
        show_anim: bool ; whether to display an animation of the simulation
        show_nth: int ; show animation every show_nth steps
        case: str ; label for the simulation case (used in plot titles)
        log: bool ; if True, record average fractions of S, I, and R per sweep
        fractions: bool, optional ; if True, write S/I/R counts to file 
                                   (“data_fractions.txt”), default = False

    returns:
        arr: np.ndarray ; final N x N lattice configuration after simulation
        average_infected: np.ndarray ; average fraction of infected sites per sweep
        average_susceptible: np.ndarray ; average fraction of susceptible sites per sweep
        average_recovered: np.ndarray ; average fraction of recovered sites per sweep

    description:
        Simulates the stochastic SIRS epidemic model on an N x N lattice over 
        `nsweeps` sweeps. At each Monte Carlo step, a random site is chosen and 
        updated according to transition probabilities:
            - Susceptible (0) → Infected (1) with probability p1 
              if at least one neighboring site is infected.
            - Infected (1) → Recovered (2) with probability p2.
            - Recovered (2) → Susceptible (0) with probability p3.
            - Immune (4) sites remain unchanged.
        Optionally animates the lattice evolution, logs average fractions per 
        sweep, and saves S/I/R counts to a text file if `fractions=True`.
    """

    if fractions:
        fractions_filename = 'data_fractions.txt'
        outfile_fractions = open(fractions_filename, 'w')
        outfile_fractions.write("<I> <S> <R>\n")
    
    
    nsteps = nsweeps*N*N

    average_infected = np.zeros(nsweeps)
    average_susceptible = np.zeros(nsweeps)
    average_recovered = np.zeros(nsweeps)
    
    
    if show_anim:
        
        if immunity_fraction == 0: 
            cmap = ListedColormap(["white", "darkred", "grey"])
        else:
            cmap = ListedColormap(["white", "darkred", "grey", "b"])
        fig = plt.figure()
        im=plt.imshow(arr, animated=True)

    


    for step in range(nsteps):
        
        #select random site
        i, j = np.random.randint(0, N, 2)
        
                
        
        # Susceptible to Infected
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

        

        # Infected to Recovered
        elif arr[i, j] == 1:
            if np.random.rand() < p2:
                arr[i, j] = 2
        
        
        
        # #immunity case
        elif arr[i, j] == 4:
             pass
        
        
        # Recovered to Susceptible
        else:
            if np.random.rand() < p3:
                arr[i, j] = 0
    

        if log and step % (N*N) == 0:
            average_infected[step//(N*N)] = np.sum(arr == 1)/2500
            average_susceptible[step//(N*N)] = np.sum(arr == 0)/2500
            average_recovered[step//(N*N)] = np.sum(arr == 2)/2500
        
        
        if step % show_nth == 0:
            

            if show_anim:
                plt.cla()
                im=plt.imshow(arr, animated=True, cmap=cmap)
                plt.title(f"SIRS Model, {case} case\nsweep: {step/(N*N)}, p1: {p1}, p2: {p2}, p3: {p3}: IF: {immunity_fraction}")
                plt.draw()
                plt.pause(0.001)
        
        if step%250 == 0 and fractions:
            outfile_fractions.write(f"{np.sum(arr == 1)} {np.sum(arr == 0)} {np.sum(arr == 2)}\n")
            

    
    
    return arr, average_infected, average_susceptible, average_recovered
















def phase_plane_plot():
    """
    args:
        None
    returns:
        None
    
    description:
        Reads simulation output data from the file "p1_p3_phase.txt" and 
        generates a heatmap of the average infected fraction <I>/N across 
        the (p1, p3) parameter plane.
    
        The input file is expected to contain columns of:
            p1, p2, p3, <I>, <I²>
        for a grid of (p1, p3) parameter combinations.
    
        The function:
            - Parses the data into arrays of p1, p3, and <I> values.
            - Maps <I> values onto a 2D grid corresponding to p1–p3 pairs.
            - Uses seaborn.heatmap() to visualise the phase plane, showing 
              how the mean infected fraction varies with p1 and p3.
            - Labels axes with p1 and p3 values and inverts the y-axis 
              for conventional plotting orientation.
    
        Produces a heatmap figure titled "<I>/N for p1 p3 plane".
    """

    filename = "p1_p3_phase.txt"
    
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
        
        mapped[i[0][0], j[0][0]] = dp[3]
    
    
    labels = np.round(all_ps, decimals=2)
    
    plt.figure()
    plt.title('<I>/N for p1 p3 plane')
    sns.heatmap(mapped.T, xticklabels=labels, yticklabels=labels)
    plt.gca().invert_yaxis()
    plt.xlabel('p1')
    plt.ylabel('p3')
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def variance_plane_plot():
    """
    args:
        None
    
    returns:
        None
    
    description:
        Reads simulation output data from "p1_p3_phase.txt" and produces 
        a heatmap of the variance in infected fraction (I) across the 
        (p1, p3) parameter plane.
    
        The input file should contain columns of:
            p1, p2, p3, <I>, <I²>
        for different combinations of infection (p1) and recovery (p3) 
        probabilities.
    
        The function:
            - Loads and parses data from the text file.
            - Calculates the variance of infection as 
              (1/2500) * (<I²> − <I>²) for each (p1, p3) pair.
            - Maps these values to a 2D grid corresponding to the 
              parameter combinations.
            - Uses seaborn.heatmap() to visualise how infection variance 
              changes across the (p1, p3) plane.
            - Inverts the y-axis for standard heatmap orientation and 
              labels the axes with p1 and p3.
    
        Produces a figure titled "variance in I for p1 p3 plane".
    """

    
    filename = "p1_p3_phase.txt"
    
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
    sns.heatmap(mapped.T, xticklabels=labels, yticklabels=labels)
    plt.title('variance in I for p1 p3 plane')
    plt.gca().invert_yaxis()
    plt.show()
    
    
    
    
    
    
    
    
    

def plot_cut():
    """
    args:
        None
    
    returns:
        None
    
    description:
        Reads simulation data from "data_from_cut.txt" and plots the 
        variance in infected fraction (I) as a function of p1 for a 
        fixed set of p2 and p3 values.
    
        The input file is expected to contain columns of:
            p1, p2, p3, <I>, <I²>, ... , error
        where the 8th column (index 7) contains precomputed jackknife 
        errors or uncertainty estimates.
    
        The function:
            - Loads the data, ignoring commented header lines (starting with 'p').
            - Computes the variance for each p1 value as 
              (1/2500) * (<I²> − <I>²).
            - Plots variance versus p1 using error bars to represent 
              uncertainty in the variance estimate.
            - Labels axes and displays the plot titled 
              "variance plot for fixed p2 p3".
    
        Produces an error-bar plot showing how infection variance changes 
        with infection probability p1 along a fixed p2–p3 cut in parameter space.
    """

    f_name = 'data_from_cut.txt'
    data = np.loadtxt(f_name, comments = 'p', dtype=float)

    p1s = data[:, 0]

    I = data[:, 3]
    Isqr = data[:, 4]
    err = data[:, 7]
    
    vals = np.zeros(len(I))
    
    for i in range(len(vals)):
        vals[i] = (1/2500)*(Isqr[i] - I[i]**2)
    
    
    plt.figure()
    #plt.plot(p1s, vals)
    plt.errorbar(p1s, vals, yerr = err, ecolor = 'r')
    plt.title('variance plot for fixed p2 p3')
    plt.xlabel('p1')
    plt.ylabel('variance')
    plt.show()





def fractions_plots():
    """
    args:
        None
    
    returns:
        None
    
    description:
        Reads time-series data of state fractions (S, I, R) from 
        "data_fractions.txt" and plots their evolution over simulation steps.
    
        The input file is expected to contain columns of:
            <I> <S> <R>
        recorded at regular intervals during the simulation.
        Each value represents the count of sites in a given state, 
        which are normalised by the total number of sites (2500) to 
        obtain fractional populations.
    
        The function:
            - Loads and normalises the data to convert counts to fractions.
            - Extracts time-series arrays for infected (I), susceptible (S), 
              and recovered (R) states.
            - Creates a scatter plot showing S, I, and R fractions versus 
              simulation step index.
            - Uses colour coding: red (I), black (S), grey (R), and labels 
              each series directly on the plot.
    
        Produces a figure titled "fraction of states" showing the 
        temporal evolution of S, I, and R fractions in the population.
    """

   
    f_name = 'data_fractions.txt'
    data = np.loadtxt(f_name, comments = '<', dtype=float, ndmin=2)
    data /= 2500
    Is = data[:, 0]
    Ss = data[:, 1]
    Rs = data[:, 2]

    x = np.array(list(range(len(Is))))
    
    plt.figure()
    plt.scatter(x, Is, c='r', s=3)
    plt.scatter(x, Ss, c='black', s=3)
    plt.scatter(x, Rs, c='grey', s=3)
    plt.xlabel('250 Step Increments')
    plt.ylabel('Fraction')
    
    plt.title("Fraction of States")
    plt.text(x[-1]*1.07, 0.55,"S", c = 'black')
    plt.text(x[-1]*1.07, 0.5,"I", c = 'r')
    plt.text(x[-1]*1.07, 0.45,"R", c = 'grey')
    plt.show()
    
    
    
    








def immunity_calculations():
    """
    args:
        None
    
    returns:
        None
    
    description:
        Runs a series of SIRS simulations across a range of immunity 
        fractions and records the resulting average infected fraction <I> 
        for each level of immunity.
    
        The function:
            - Defines a set of immunity fractions (`fractions`) ranging 
              from 0 to 1, with finer spacing at low immunity.
            - For each fraction f:
                • Initialises an N×N lattice with fraction f of immune sites.
                • Runs an equilibration phase (no data recorded).
                • Runs a measurement phase and records the average infected 
                  fraction <I> over time.
            - Writes results to "immunity_data.txt" in the format:
                  I_fi   <I>
              where I_fi is the fraction of immune sites and <I> is the mean 
              infected fraction after equilibration.
    
        The simulation uses fixed parameters:
            N = 50, p1 = p2 = p3 = 0.5
            100 equilibration sweeps and 100 measurement sweeps.
    
        Produces a data file showing how increasing immunity affects the 
        steady-state infection level in the SIRS model.
    """

    fractions = np.concatenate((np.arange(0, 0.4, .005), np.arange(0.4, 1, 0.1)))
    print(fractions)
    N = 50
    p1 = 0.5
    p2 = 0.5
    p3 = 0.5
    
    
    show_anim = False
    
        
    show_nth = 100
    case = 'Default'
    
    equilibration_sweeps = 100
    measurement_sweeps = 100

    filename = 'immunity_data.txt'
    outfile_im = open(filename, 'w')
    outfile_im.write('I_fi <I>\n')


    for f in fractions:
        print(f)
        arr = initialise(N, f)
        
        #equilibrate
        nsweeps = equilibration_sweeps
        arr, ignore, ignore, ignore = simulate(arr, p1, p2, p3, N, nsweeps, show_anim, show_nth, case, log = False, fractions = f)
        
        
        nsweeps = measurement_sweeps
        arr, av_infected, average_susceptible, average_recovered = simulate(arr, p1, p2, p3, N, nsweeps, show_anim, show_nth, case, log = True, fractions = f)
        
        outfile_im.write(F"{f} {np.mean(av_infected)}\n")
        

        
        

    
def immunity_plots():
    """
    args:
        None
    
    returns:
        None
    
    description:
        Reads the simulation output from "immunity_data.txt" and plots 
        the relationship between the fraction of immune sites and the 
        average infected fraction <I>.
    
        The input file is expected to contain two columns:
            I_fi   <I>
        where I_fi is the fraction of immune sites and <I> is the 
        corresponding average infected fraction recorded from the simulation.
    
        The function:
            - Loads the data from the text file.
            - Extracts arrays for immunity fractions and mean infected values.
            - Creates a scatter plot of <I> versus fraction of immunity.
            - Labels the axes and titles the figure as 
              "average infected sites against fraction of vaccinations".
    
        Produces a figure that visualises how increasing immunity reduces 
        the average infection in the SIRS model.
    """

    
    data = np.loadtxt('immunity_data.txt', comments = 'I', dtype=float)
    fs = data[:, 0]
    Is = data[:, 1]
    
    plt.figure()
    plt.scatter(fs, Is, s = 5)
    plt.xlabel("fraction of immunity")
    plt.ylabel("average infected sites")
    plt.title("average infected sights against fraction of vaccinations")
    plt.show()















def main():
    """
    args:
        None (all parameters are handled via command-line arguments)
    
    returns:
        None
    
    description:
        Entry point for the SIRS lattice simulation. Handles command-line 
        arguments to control simulation behaviour, data collection, and 
        plotting.
    
        Key functionalities:
            - Runs single-case or automated simulations across a range of 
              parameters (p1, p2, p3).
            - Supports predefined cases: 'absorbing', 'dynamic_eq', 'cyclic'.
            - Allows setting lattice size, number of sweeps, and display of 
              animation.
            - Automates data collection for phase-plane exploration or cuts 
              through parameter space.
            - Records simulation outputs to text files for further analysis:
                • 'p1_p3_phase.txt' for full parameter plane
                • 'data_from_cut.txt' for cuts at fixed p2, p3
            - Supports plotting functions:
                • Phase-plane average infection (<I>/N)
                • Phase-plane variance ((<I²> - <I>²)/N)
                • Fractions of S, I, R over time
                • Variance along a cut
                • Immunity vs average infection

        The function orchestrates the simulation workflow, calls the 
        simulate() function for lattice updates, logs results, and triggers 
        plotting as requested by the user.
    """

    
    cases = {
        'absorbing' : np.array([0.1, 0.9, 0.9], dtype=np.float64),
        'dynamic_eq' : np.array([0.5, 0.5, 0.5], dtype=np.float64),
        'cyclic' : np.array([0.8, 0.1, 0.01], dtype=np.float64)
    }
    
    parser = argparse.ArgumentParser(description="SIR Model Simulation with CLI Inputs")
    parser.add_argument("-a", "--auto", action='store_true', help="full auto mode for data collection")
    parser.add_argument("-ac", "--auto_cut", action = 'store_true', help = "auto mode for data collection along specified cut")
    parser.add_argument("-N", "--N", type=int, default=50, help="Size of the lattice")
    parser.add_argument("-s", "--show_anim", action= "store_true", help="Show animation of the simulation")
    parser.add_argument("-ns", "--nsweeps", type=int, default=500, help="Number of sweeps")
    parser.add_argument("-p1", "--p1", type=float, default=0.5, help="Probability of S -> I transition")
    parser.add_argument("-p2", "--p2", type=float, default = 0.5, help="Probability of I -> R transition")
    parser.add_argument("-p3", "--p3", type=float, default = 0.5, help="Probability of R -> S transition")
    parser.add_argument("-c", "--case", choices=["absorbing", "dynamic_eq", "cyclic"], help="Choose a predefined case")
    parser.add_argument("-nth", "--show_nth", type=int, default=2500, help="Show every nth step")
    parser.add_argument("-phase", "--p1_p3_phase", action = 'store_true', help = "Plots <I> / N phase across p1, p3 plane")
    parser.add_argument("-var", "--phase_variance", action = 'store_true', help = "Plots ( <I2> - <I>2 ) / N phase across p1, p3 plane")
    parser.add_argument("-f", "--fractions", action = 'store_true', help = "Plots the fractions of each state across the simulation")
    parser.add_argument("-pc", "--plot_cut", action = 'store_true', help = "Plots the variance of cut simulation")
    parser.add_argument("-if", "--immunity_fraction", type = float, default = 0, help = "fraction of immune/vaccinated agents")
    parser.add_argument("-ai", "--auto_immunity", action = 'store_true', help = "automates imunity against fractions plot")
    parser.add_argument("-pi", "--plot_immunity", action = 'store_true', help = "Plots data from immunity_data.txt")


    args = parser.parse_args()
    
    global immunity_fraction
    immunity_fraction = args.immunity_fraction
    
    
    if args.plot_immunity:
        immunity_plots()
        return
    
    if args.auto_immunity:
        immunity_calculations()
        return

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
        measurement_sweeps = 10000
        
        #unused
        show_nth = 100
        case = 'Default'
        
        cut_filename = 'data_from_cut.txt'
        outfile_cut = open(cut_filename, 'w')
        outfile_cut.write("p1 p2 p3 <I> <I^2> <S> <R> <I_err>\n")
        
        

    elif args.auto:
        
        # all possible ps
        all_ps = np.arange(0, 1.05, 0.05)
        p2s = np.array([0.5], dtype=float)
        p1s = all_ps
        p3s = all_ps
        
        N = 50
        show_anim = args.show_anim
        equilibration_sweeps = 100
        measurement_sweeps = 1000
        
        #unused
        show_nth = 100
        case = 'Default'
  
        filename = "p1_p3_phase.txt"
        outfile_plane = open(filename, 'w')
        outfile_plane.write("p1 p2 p3 <I> <I^2> <S> <R> <I_err>\n")
    
    
    
    

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
            
            
                arr = initialise(N, args.immunity_fraction)
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
                    outfile_plane.write(f"{p1} 0.5 {p3} {np.mean(av_infected)} {np.mean(av_infected**2)} {np.mean(average_susceptible)} {np.mean(average_recovered)} {jackknife(av_infected)}\n")
                    print(f"{p1} 0.5 {p3} {np.mean(av_infected)} {np.mean(av_infected**2)} {np.mean(average_susceptible)} {np.mean(average_recovered)}\n")
                
                if args.auto_cut:
                    outfile_cut.write(f"{p1} 0.5 {p3} {np.mean(av_infected)} {np.mean(av_infected**2)} {np.mean(average_susceptible)} {np.mean(average_recovered)} {jackknife(av_infected)}\n")
                    print(f"{p1} 0.5 {p3} {np.mean(av_infected)} {np.mean(av_infected**2)} {np.mean(average_susceptible)} {np.mean(average_recovered)}\n")


    
    if args.fractions:
        fractions_plots()
main()


#if elif in simulate for 
