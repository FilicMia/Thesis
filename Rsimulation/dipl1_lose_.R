#NE VALJA,bolje kada se modelira odvojeno.

library("MASS","matrixcalc","base")
library("Matrix", lib.loc="C:/Program Files/R/R-3.3.2/library")
#pseudo-udaljenosti
c = 299792458
R = read.csv('pseudoranges.txt', header = TRUE);
R <- as.matrix(R$R)
#učitaj koordinate satelita
S = read.csv('satellites.txt')
S <- as.matrix(S)
x_0 = c(1,1,1,1)
delt = c(3,3,3,3)

unutar = append(S,R)
RS = matrix(unutar,4,4) # [x,y,z,c*dT]

while(norm(t(delt)) > 1){
  xRadni = x_0
  xRadni[4] = c*xRadni[4] #zadnji s c
  
  P = t(apply(RS, 1, function(x) (x-xRadni))) #x-x_0
  PP = P
  PP[,4] = - c*P[,4] #4. red s -c
  
  PX = P*P
  PX = P%*%c(-1,-1,-1,1)
  
  delt <- -(1/2)*svd.inverse(PP)%*%PX
  x_0 = x_0 + delt #(x,y,z,c*dT)
  print(delt)  
}