Checkpoint 03 - Modelling and Visualisation

Cahn-hilliard equation solved in file 'ch.py', 
    further details in argparser if needed

BVP Poissons equation solved in file 'poisson.py', 
    only using jacobi method further details in argparser if needed

BVP Magnet - solved in 'cb.py' 
    uses gauss-seidel (SOR only) and jacobi algorithms, use -l argument if user is wanting to
    see plots of newly run simulation (nb this overwrites previous runs)
    further details provided in argparser

    for purely the gauss-seidel method - see 'magnet.py', note this 
    doesnt provide vector plots for this method




Extra : the file 'testmag.py' solves the magnet problem with the SOR, gauss siedel, and
    jacobi algorithms, the file is messy but it provides all plots including vector plots
    and functions of distance plots which the aforementioned files lack. nb ensure to use -l
    argument here too to plot the data of the current run. (essentially combines cb.py and magnet.py)


for all plots that take too long to run (particularly with GS) .pngs are labelled
and stored in the file.