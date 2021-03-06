rm(list=ls())

library(dplyr)
library(MASS)

npass <- 0 # po�etni postav brojila iteracija

x <- read.csv('satellites.txt', header = FALSE) # u�itavanje polo�aja satelita, stupci su (x,y,z)
# x <- as.matrix(x)
pr <- read.csv('pseudoranges.txt', header = FALSE) # u�itavanje izmjerenih pseudoudaljenosti
pr <- as.vector(pr[,1])
pr <- t(t(pr))

dx <- c(10,10,10, 10) # postav po�etnih uvjeta za iteraciju - razlika susjednih iteracija (zaustavlja iteraciju)
X <- t(t(c(13, 0, 44, 5))) # postav po�etnih uvjeta za iteraciju (procjena polo�aja)

# Definicija matrica i fizkalnih konstanti
G <- matrix(nrow = length(t(pr)), ncol = length(x)+1)
R <- c(1,1)
C <- diag(length(t(pr))) # matrica kovarijanci za Weighted LS solution - u po�etku postavljena kao jedini�na matrica
                         # za regular position estimation (pretpostavljena potpuna kompenzacija pogre�aka)
c <- 2.99792458E+08 # brzina svjetlosti [m/s], po GPS standardu
eps <- 3.0 # najve�a prihvatljiva pogre�ka komponente odre�ivanja polo�aja [m], eps = max(eps(x), eps(y), eps(z))

while(eps < max(abs(dx[1]), abs(dx[2]), abs(dx[3]), abs(dx[4]))){ # iteracija - sve dok sve pogre�ke po komponentama ne budu manje od eps
for(i in 1:length(t(pr))){
  R[i] <- sqrt((x[i,1] - X[1])^2 + (x[i,2] - X[2])^2 + (x[i,3] - X[3])^2)
  pr[i] <- t(t(pr[i])) - R[i] #pr = rho
  
  for(j in 1:(length(x))){
    G[i,j] <- (X[j] - x[i,j])/R[i]
  }
  
  G[i, length(x)+1] <- c
  
  }
npass <- npass + 1
print(npass)
print(max(abs(dx[1]), abs(dx[2]), abs(dx[3]), abs(dx[4])))
	print(G)

# Postupak procjene polo�aja
ginv = solve(G)
dx <- ginv %*%  C %*% pr # !! PROBLEM!!! singularna matrica dobivena sa solve??
  
# W <- t(G) %*% G # prvo mogu�e alternativno rje�enje
# W <- qr(W) 
# dx <- W$qr %*% t(G) %*% C %*% pr 
  
# dx <- solve(G, tol = 1E-50) %*% pr # drugo mogu�e alternativno rje�enje 
X <- X + dx
}

#TO DO: data collection for graphical presentation of:
# - efficiency of the iteration algorithm
# - position estimation accuracy
# - interpretation
