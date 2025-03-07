import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d





def initial_condition(N):
    arr = np.zeros((N, N))
    
    noise = np.random.rand(N, N)*0.1
    arr += noise
    
                
    return arr



def testing():

    
    N = 50
    a = 0.1
    M = 0.1
    kappa = 0.1

    dt = 1
    dx = 1
    
    kernel = np.array([[0, 1, 0],
                       [1, 0, 1],
                       [0, 1, 0]])
    
    
    
    phi = initial_condition(N)
    

    conv_phi= convolve2d(phi, kernel, mode='same', boundary='wrap')
    mu = -a*phi + a*phi**3 - kappa*conv_phi - 4*phi
    
    conv_mu = convolve2d(mu, kernel, mode='same', boundary='wrap')
    phi_new = phi + ((M*dt)/(dx**2)) * (conv_mu - 4*mu)
    
    print(phi_new)
    
    








def test_loop():
    
    
        
    N = 100
    a = 0.1
    M = 0.1
    kappa = 0.1

    dt = 1
    dx = 1
    
    kernel = np.array([[0, 1, 0],
                        [1, 0, 1],
                        [0, 1, 0]])
    

    
    phi = initial_condition(N)
    step = 0
    while True:
        step +=1

        conv_phi= convolve2d(phi, kernel, mode='same', boundary='wrap')
        mu = -a*phi + a*phi**3 - kappa*(conv_phi - 4*phi)/(dx**2)
        
        conv_mu = convolve2d(mu, kernel, mode='same', boundary='wrap')
        phi +=  ((M*dt)/(dx**2)) * (conv_mu - 4*mu)
                
        if step % 100== 0:
            plt.cla()
            plt.imshow(phi, animated=True)
            plt.draw()
            plt.pause(0.01)
            
            
            print(np.sum(phi))
            
        
            

test_loop()
    
    
