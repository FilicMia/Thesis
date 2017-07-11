#dipl1 isprobati na ovim podatcima isto
library("MASS","matrixcalc","base", "Matrix")
#pseudo-udaljenosti
c <- 2.99792458E+08 # brzina svjetlosti [m/s], po GPS standardu
R = read.csv('pseudoranges5a.txt', header = FALSE);
R <- as.matrix(R[,1])
#učitaj koordinate satelita
S = read.csv('satellites5.txt', header = FALSE)
S <- as.matrix(S)
x_0 = c(1,1,1,1)
delt = c(3,3,3,3)

unutar = append(S,R)
RS = matrix(unutar,dim(S)[1],dim(S)[2]+1) # [x,y,z,c*dT]

while(norm(t(delt)) > 1){
  
  P = t(apply(RS, 1, function(x) (x-x_0))) #x-x_0
  PP = P
  PP[,dim(S)[2]+1] = - c*P[,dim(S)[2]+1] #4. red s -c
  
  A = -2*PP
  
  #QR od A!!!!
  x_0 = x_0 + delt #(x,y,z,c*dT)
  print(delt)  
}