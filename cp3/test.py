import numpy as np

N = 5

arr = np.zeros((N, N, N))

for i in range(N):
    for j in range(N):
        for k in range(N):
            
            
            if (i+j+k) % 2 == 0:
                
                arr[i, j, k] = 1
                    

# print(arr)


evens = np.empty((0, 3), dtype=int)  
odds = np.empty((0, 3), dtype=int)  

for i in range(N):
    for j in range(N):
        for k in range(N):
            
            if (i+j+k) % 2 == 0:
                evens = np.vstack([evens, [i, j, k]])
            
            else:
                odds = np.vstack([odds, [i, j, k]])
            
            
print(evens)