import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp




def initialise(N):
    matrix = np.random.choice([True, False], (N,N))
    print(matrix)


def simulate():
    
    #Params
    N=50
    
    
    matrix = initialise(N)
    