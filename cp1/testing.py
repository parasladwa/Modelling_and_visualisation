# import sys
# import math
# import random
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# import time



# def set_up(N):
#     #randomly set up array
#     arr = np.zeros([N, N], dtype=int)
    
#     for i in range(0, len(arr)):
#         for j in range(0, len(arr[0])):
#             arr[i][j] = random.choice([-1, 1])
    
#     return arr


# arr = set_up(50)
# #copied
# fig = plt.figure()
# im=plt.imshow(arr, animated=True)

# for i in range(1000):
#     arr = set_up(50)
  
#     plt.cla()
#     im=plt.imshow(arr, animated=True)
#     plt.draw()
#     plt.pause(0.0001)





# from numba import jit, cuda 
# import numpy as np 
# # to measure exec time 
# from timeit import default_timer as timer    
  
# # normal function to run on cpu 
# def func(a):                                 
#     for i in range(10000000): 
#         a[i]+= 1      
  
# # function optimized to run on gpu  
# @jit(target_backend='cuda')                          
# def func2(a): 
#     for i in range(10000000): 
#         a[i]+= 1
# if __name__=="__main__": 
#     n = 10000000                            
#     a = np.ones(n, dtype = np.float64) 
      
#     start = timer() 
#     func(a) 
#     print("without GPU:", timer()-start)     
      
#     start = timer() 
#     func2(a) 
#     print("with GPU:", timer()-start) 








for i in range(10):
    if i == 5:
        print('here')
        break #continue
        
        
    print(i)







    
    
    
    
    
    
    