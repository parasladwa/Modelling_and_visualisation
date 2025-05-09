
HOW TO RUN

    all code is within exam.py,
    using argparse each part can be run independantly in main():
    
        a-  simulation

        b- no code to run as it is an image based question

        c- runs full simulations and plots variance against initial phi_0
            the image is already saved in the file

        d- separate simulation function, still can be called within main
            simulation includes advection term.






part b discussion
    for initial phi0 == 0 (part i) we can see there is an

    for initial phi0 == -0.1 (part ii) we observe nucleation favouring one chemical

    for initial phi0 == 0.1 (part iii) we observe nucleation favouring the opposite
                    chemical as in part ii
    
    for all cases this resembles the cahn hilliard equation


part c discussion - refer to images and datafile as well as below 
    
    the variance allows us to see how strongly the separation of the chemicals are
    we can see that this begins high and drops significantly, the maximal gradiant of this
    drop being at initial phi0 == 0.1

    for smaller initial phi0 (approx 0) we have much stronger boundaries relative to larger
    initial phi0 (phi0 approx 0.25) this is evident in the plot of variance too

    





part d - refer to images as well as below

    for the smaller two advection terms (v_0 = 0.001, 0.01, images i and ii)
    we can see negligable motion

    for the largest case, v_0 == 0.1, there is prominant motion upon two channels of
    the phi lattice, moving in opposite directions, this is more visible during the
    realtime simulation in motion as opposed to the image.


