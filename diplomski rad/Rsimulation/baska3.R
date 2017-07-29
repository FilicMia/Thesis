#source: http://math.tut.fi/courses/MAT-45806/mathematics_and_methods_for_positioning_2008.pdf
#The orthogonal projection matrix P corresponding to matrix x is defined as P=x(xTx)
library(LDRTools)
#uključujemo BLUE

#učitaj podatke
c = 299792458
R = read.csv('pseudoranges.txt', header = TRUE);
R <- as.matrix(R$R)
#read saletites coordinates
S = read.csv('satellites.txt')
S <- as.matrix(S)
x_0 = c(0,0,0)
error = c(3,3,3)
while(norm(t(error)) > 2){
  preRhoVector = t(apply(S, 1, function(x) (x-x_0))) #x-x_0
  rho = t(t(apply(preRhoVector,1, function(x) norm(t(x),'F'))))
  A = t(apply(preRhoVector,1, function(x) c(x/norm(t(x),'F'), 1) ))
  b = R-rho
  
  H = B2P(A)
  M = diag(length(R))-H
  #procjenitelj za sigma = E(error)
  sigma = 1/length(R)*(t(b)%*%M%*%b) #primjenjena statistika
  
  #http://math.tut.fi/courses/MAT-45806/mathematics_and_methods_for_positioning_2008.pdf#27
  C = solve(sigma)[1,1] * diag(length(R))
  #BLUE - http://www.rtklib.com/prog/manual_2.4.2.pdf#158
  AT = t(A)
  print(AT)
  ATACATC = solve(C%*%A)%*%solve(AT)%*%AT%*%C 
  
  error <- ATACATC%*%b
  x_0 = x_0 + error[0:3,]
  print(norm(error))
}


