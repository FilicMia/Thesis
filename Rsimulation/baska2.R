#Simplefied position estimation calculator.
#Let T = t+tau - universal clock time at point of the receive of the signal = true receive clock time + clock bias
#T_s = t_s + tau_s - Universal clock time at point of the transmition of the signal = true satellite time + satellite clock bias(due to the satellite)
#rho_s(t, t_s) is the range from receiver (at receive time) to the satellite (at transmit time)
#R = rho_s(t,t_s) + c(tau-tau_s) (= (T-T_s)c )
#
#This  assumes  that  the  inverse  to 
#(A^T)A 
#exists.  For  example, m≥4is  a  necessary  (but  not sufficient)  condition.    Problems  can  exist  if,  for  example,  a  pair  of  satellites  lie  in  the  same 
#line  of  sight,  or  if  the  satellites  are  all  in  the  same  orbital  plane.  In  almost  all  practical 
#situations, m≥5is  sufficient.    Alternatively,  one  parameter  could  be  left  unestimated  (e.g.,                                                                    the height could be fixed to sea-level for a boat)
#

#read R (pseudo ranges calculated from PRN code) - observed
R = read.csv('pseudoranges.txt', header = TRUE);
R <- as.matrix(R$R)
#read saletites coordinates
S = read.csv('satellites.txt')
S <- as.matrix(S)
#linearisation.. b = Ax + error
# solution due to the least square solution of x is 
# x = ((A#t)A)^(-1)A#t * b
# where A = colums of partial derivations od P, each row, one satellite
#b = R-rho
x_0 = c(8,099,0) # 
c = 299792458
# c(x0-x_s,y0-y_s,z0-z_s,tau0-tau_s)
#A does not change, no matter that is it x-xo/rho and xo changes as well as rho
error = c(3,3,3)

while(norm(t(error)) > 2){
  preRhoVector = t(apply(S, 1, function(x) (x-x_0))) #x-x_0
  rho = t(t(apply(preRhoVector,1, function(x) norm(t(x),'F'))))
  A = t(apply(preRhoVector,1, function(x) c(x/norm(t(x),'F'), 1) )) #x_0[4] = detla_t *c
  b = R-rho
  AT = t(A)
  ATAAT = solve(A)%*%solve(AT)%*%AT # = solve(AT%*%A)%*%AT <- MLE procjenitelj
  
  error <- ATAAT%*%b
  x_0 = x_0 + error[0:3,]
  print(norm(error))
}
#pretpostavljam da je kovarijanca prevelika pa da ne konvergira.

#The orthogonal projection matrix P corresponding to matrix x is defined as P=x(xTx)


