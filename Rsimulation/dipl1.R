library("MASS","matrixcalc")
#pseudo-udaljenosti
c = 299792458
R = read.csv('pseudoranges.txt', header = TRUE);
R <- as.matrix(R$R)
#učitaj koordinate satelita
S = read.csv('satellites.txt')
S <- as.matrix(S)
x_0 = c(0,0,0,0)
delt = c(3,3,3,3)

unutar = append(S,R)
RS = matrix(unutar,4,4) # [x,y,z,dT]

while(norm(t(delt)) > 1){
  P = t(apply(RS, 1, function(x) (x-x_0))) #x-x_0
  PP = P
  PP[,4] = - c*P[,4] #4. red s -c
  
  PX = P*P
  PX = PX%*%c(-1,-1,-1,1)
  
  delt <- -(1/2)*svd.inverse(PP)%*%PX
  x_0 = x_0 + delt #(x,y,z,dT)
  print(delt)  
}

#Konvergira, ali treba provjeriti je li rješenje dobro!!!!
# [,1]
# [1,]  4.369589e+06
# [2,] -3.346442e+06
# [3,]  3.162960e+06
# [4,]  6.254601e-04
# [,1]
# [1,]  3.208789e+05
# [2,] -8.578050e+04
# [3,]  7.142044e+05
# [4,]  9.782826e-06
# [,1]
# [1,]  5.097270e+03
# [2,] -1.362652e+03
# [3,]  1.134538e+04
# [4,]  2.585285e-09
# [,1]
# [1,]  1.287138e+00
# [2,] -3.439509e-01
# [3,]  2.864487e+00
# [4,]  1.164256e-10
# [,1]
# [1,]  4.683763e-08
# [2,] -2.306779e-08
# [3,]  1.503124e-07
# [4,]  1.164253e-10

# x_0
# [,1]
# [1,]  4.695567e+06 x
# [2,] -3.433585e+06 y
# [3,]  3.888512e+06 z
# [4,]  6.352458e-04 d_T

