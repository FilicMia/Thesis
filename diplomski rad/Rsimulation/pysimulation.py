from numpy import *

# generating a random overdetermined system
A = [
[22887118, -43707444, -18568101, -1.669527e+17],
[-17699729, -30423008, -39675983, -1.670857e+17],
[25759894, -16855821, -43459194, -1.655726e+17],
[12447714, -51004833, 7679259, -1.665416e+17],
[-20852090,
-43656310,
22151333,
-1.658202e+17],
[-3022873,
-47396299,
-23299744,
-1.677310e+17],
[-39662313,
-13121244,
-34412557,
-1.654622e+17]]
b = [0]*len(A)

x_lstsq = linalg.lstsq(A,b)[0] # computing the numpy solution

Q,R = linalg.qr(A) # qr decomposition of A
Qb = dot(Q.T,b) # computing Q^T*b (project b onto the range of A)
x_qr = linalg.solve(R,Qb) # solving R*x = Q^T*b

# comparing the solutions
print ('qr solution')
print (x_qr)
print ('lstqs solution')
print (x_lstsq)
