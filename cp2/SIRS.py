import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def initialise(N):
    
    arr = np.random.choice([0, 1, 2], size=(N, N))
    
    return arr



def main():
    # S0 I1 R2
    N = 50
    show_anim = True
    show_nth = 100
    arr = initialise(N)
    nsweeps = 100
    nsteps = nsweeps*N*N

    
    p1, p2, p3 = 0.5, 0.5, 0.5
    
    
    if show_anim:
        cmap = ListedColormap(["white", "red", "grey"])
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
        if arr[i, j] == 1:
            if np.random.rand() < p2:
                arr[i, j] = 2
        
        # R2 -> S0
        if arr[i, j] == 2:
            if np.random.rand() < p3:
                arr[i, j] = 0
        

        
        if show_anim and step % (show_nth%(N*N)) == 0:
            plt.cla()
            im=plt.imshow(arr, animated=True, cmap=cmap)
            plt.text(1, -1, f"{step}")
            plt.draw()
            plt.pause(0.01)
main()










