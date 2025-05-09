import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve
kernel2d = np.array([[0,1,0],
                   [1,0,1],
                   [0,1,0]])


def poisson_solver_2d(phi,mu,dx,dt,factor):
    new_phi = np.zeros_like(phi)
    new_phi += phi + (factor * dt / (dx**2)) * (convolve(mu, kernel2d, mode='wrap')-4*mu)
    return new_phi

def chemical_solver(a,b,c,q,p,D,dt):
    term1 = poisson_solver_2d(a,a,1,dt,D)
    term2 = q*a*(1-(a+b+c))-p*a*c
    return term1 + term2


def update_function(chemical_a,chemical_b,chemical_c,chemical_triple,typefield,q,p,D,dt):
    new_chemical_a = chemical_a.copy()
    new_chemical_b = chemical_b.copy()
    new_chemical_c = chemical_c.copy()
    # update a
    new_chemical_a = chemical_solver(chemical_a,chemical_b,chemical_c,q,p,D,dt)
    # update b
    new_chemical_b = chemical_solver(chemical_b,chemical_c,chemical_a,q,p,D,dt)
    # update c
    new_chemical_c = chemical_solver(chemical_c,chemical_a,chemical_b,q,p,D,dt)
    # update triple
    chemical_triple = 1-(chemical_a+chemical_b+chemical_c)
    # update typefield
    typefield = np.argmax(np.stack([chemical_a, chemical_b, chemical_c, chemical_triple], axis=0), axis=0)
    return new_chemical_a, new_chemical_b, new_chemical_c, chemical_triple, typefield



def main():

    gsize = 50
    q = 1
    p = 1
    D = .5
    k = 1/3

    dt = .1


    chemical_a = np.random.uniform(0, k, (gsize, gsize))
    chemical_b = np.random.uniform(0, k, (gsize, gsize))   
    chemical_c = np.random.uniform(0, k, (gsize, gsize))

    chemical_triple = np.zeros((gsize,gsize))
    typefield = np.zeros((gsize,gsize))
    chemical_a, chemical_b, chemical_c, chemical_triple, typefield = update_function(chemical_a, chemical_b, chemical_c, chemical_triple, typefield, q,p,D,dt)

    
    print("begin")
    colors = ['blue', 'green', 'red', 'yellow']
    cmap = plt.matplotlib.colors.ListedColormap(colors)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = plt.matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    
    # Create a legend for the key
    labels = ['Chem A', 'Chem B', 'Chem C', 'Triple']
    patches = [plt.matplotlib.patches.Patch(color=colors[i], label=labels[i]) for i in range(len(colors))]
    
    while True:
        chemical_a, chemical_b, chemical_c, chemical_triple, typefield = update_function(chemical_a, chemical_b, chemical_c, chemical_triple, typefield, q, p, D, dt)
        plt.cla()
        #print("\n \n chemical_a", np.sum(chemical_a),"\n chemical_a", np.sum(chemical_b),"\n chemical_c", np.sum(chemical_c),"\n Triple", np.sum(chemical_a))
        
        im = plt.imshow(typefield, cmap=cmap, norm=norm, animated=True)
        plt.legend(handles=patches, loc='upper right', title="Key")
        plt.draw()
        plt.pause(0.01)

	
main()


def main2():

    gsize = 50
    q = 1
    p = 0.5
    D = 1
    k = 1/3

    dt = .1


    chemical_a = np.random.uniform(0, k, (gsize, gsize))
    chemical_b = np.random.uniform(0, k, (gsize, gsize))   
    chemical_c = np.random.uniform(0, k, (gsize, gsize))

    chemical_triple = np.zeros((gsize,gsize))
    typefield = np.zeros((gsize,gsize))
    chemical_a, chemical_b, chemical_c, chemical_triple, typefield = update_function(chemical_a, chemical_b, chemical_c, chemical_triple, typefield, q,p,D,dt)

    
    print("begin")
    colors = ['blue', 'green', 'red', 'yellow']
    cmap = plt.matplotlib.colors.ListedColormap(colors)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = plt.matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    
    # Create a legend for the key
    labels = ['Chem A', 'Chem B', 'Chem C', 'Triple']
    patches = [plt.matplotlib.patches.Patch(color=colors[i], label=labels[i]) for i in range(len(colors))]
    
    while True:
        chemical_a, chemical_b, chemical_c, chemical_triple, typefield = update_function(chemical_a, chemical_b, chemical_c, chemical_triple, typefield, q, p, D, dt)
        plt.cla()
        print("\n \n chemical_a", np.sum(chemical_a),"\n chemical_a", np.sum(chemical_b),"\n chemical_c", np.sum(chemical_c),"\n Triple", np.sum(chemical_a))
        
        im = plt.imshow(typefield, cmap=cmap, norm=norm, animated=True)
        plt.legend(handles=patches, loc='upper right', title="Key")
        plt.draw()
        plt.pause(0.01)

	
#main2()