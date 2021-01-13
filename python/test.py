import numpy as np
A = np.zeros((2,5))
for i in range(5):
    for j in range(2):
        A[j,i] = i**(j+1)
print(A)
print(A.sum())
print(A.size)
print(A.sum()/A.size)
