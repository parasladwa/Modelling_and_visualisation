import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import time


#dummy
type_field = None
N = None
step = None
shownth = None
cmap = None
a = None
D = None
p = None
b = None
c = None
dt = None
abc = None
q = None
dx = None
phi = None
kappa = None




# show steps animation cmap also

from matplotlib.colors import ListedColormap
cmap = ListedColormap(['gray', 'red', 'green', 'blue'])

if step%shownth == 0: 
    plt.cla()
    plt.imshow(type_field, animated=True, cmap = cmap)
    plt.text(0, -.1*N, f"Step: {step}")
    plt.draw()
    plt.pause(0.01)
    
    
    
    
    
    
#laplacian convolve convolve2d 
# in the form dmatrix / dt = something laplacian A


from scipy.signal import convolve2d

kernel = [[0, 1, 0],
            [1, 0, 1],
            [0, 1, 0]]

conv_a = convolve2d(a, kernel, mode = 'same', boundary = 'wrap')
a += (D*(conv_a-4*a)/(dx**2) + q*a*(abc) - p*a*c)*dt

conv_phi= convolve2d(phi, kernel, mode='same', boundary='wrap')
mu = -a*phi + a*phi**3 - kappa*(conv_phi - 4*phi)/(dx**2)



#argparse arguments args
import argparse
parser = argparse.ArgumentParser(description="Solve the Cahn-Hilliard equation")

parser.add_argument("-a", "--a", type=float, default=0.1, help="Parameter a")
parser.add_argument("-nth", "--show_nth", type=int, default=100, help="Show every nth step")

parser.add_argument("-s", "--show_anim", action="store_true", help="show animation")
parser.add_argument("-c", "--case", type=str, default = 'zeros', help = "choose a case from {zeros, halfs}")
args = parser.parse_args()

N = args.N
a = args.a
