import numpy as np
import matplotlib.pyplot as plt

def initialise(N):
    
    arr = np.random.choice([0, 1, 2], size=(N, N))
    
    return arr



def main():
    
    N = 50
    show_anim = True
    arr = initialise(N)
    nsweeps = 10
    nsteps = nsweeps*N*N
    
    p1, p2, p3 = 0.5, 0.5, 0.5
    
    
    if show_anim:
        fig = plt.figure()
        im=plt.imshow(arr, animated=True)
        
    
    for step in range(nsteps):
        
        i, j = np.random.randint(0, N, 2)
        
        print(i, j)
main()










