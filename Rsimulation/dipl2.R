#korišteni novi podatci
library("MASS","matrixcalc")
#pseudo-udaljenosti
c <- 2.99792458E+08 # brzina svjetlosti [m/s], po GPS standardu
R = read.csv('pseudoranges5a.txt', header = FALSE);
R <- as.matrix(R[,1])
#učitaj koordinate satelita
S = read.csv('satellites5.txt', header = FALSE)
S <- as.matrix(S)
x_0 = c(1,1,1,1) #[x,y,z,d_T] d_t se kasnije množi sa c da bi se oduzeo od [x_i,y_i,z_i,d] 
delt = c(3,3,3,3)

realPosition = c(918074.1038,5703773.539,2693918.9285,0)

unutar = append(S,R)
RS = matrix(unutar,dim(S)[1],dim(S)[2]+1) # [x_i,y_i,z_i,d] 

iter = 0
niter = 100
err <- c(11,11,11,11)
#while(norm(t(delt)) > 1){
start.time <- Sys.time()

while(iter < niter){
  x_ = c(x_0[1:3],x_0[4]*c)
  P = t(apply(RS, 1, function(x) (x-x_))) #x-x_0
  PP = P
  PP[,dim(S)[2]+1] = - c*P[,dim(S)[2]+1] #4. red s -c
  
  PX = P*P
  PX = PX%*%c(-1,-1,-1,1)
  
  delt <- -(1/2)*svd.inverse(PP)%*%PX
  x_0 = x_0 + delt #(x,y,z,c*dT)
  
  cat(c(iter, delt[1:3]),' \r',file="razmakIteracija.txt", append=TRUE) # upisivanje vrijednosti dx radi kasnije analize brzine i to?nosti postupka
  err <- x_0 - realPosition
  cat(c(iter, err[1:3]),' \r',file="stvarnoOdstupanje.txt", append=TRUE)
  
  iter = iter +1
  print(delt)  
}

end.time <- Sys.time()

timediff <- end.time - start.time

d_iter <- read.csv('razmakIteracija.txt', header = FALSE, sep = '')
err <- read.csv('stvarnoOdstupanje.txt', header = FALSE, sep = '')

iter <- d_iter$V1
d.x <- d_iter$V2
d.y <- d_iter$V3
d.z <- d_iter$V4

plot(iter, log(abs(d.x)), type = 'l', col = 'red', main = c('LSA Time of execution in [s]=', round(timediff, digits = 2)), xlab = 'No. of iterations', ylab = 'log of d_X components')
lines(iter, log(abs(d.y)), type = 'l', col = 'green')
lines(iter, log(abs(d.z)), type = 'l', col = 'blue')

plot(iter, (abs(d.x)), type = 'l', col = 'red', main = c('LSA Time of execution in [s]=', round(timediff, digits = 2)), xlab = 'No. of iterations', ylab = 'd_X components')
lines(iter, (abs(d.y)), type = 'l', col = 'green')
lines(iter, (abs(d.z)), type = 'l', col = 'blue')

iter <- err$V1
xx <- err$V2
yy <- err$V3
zz <- err$V4

plot(iter, log(abs(xx)), type = 'l', col = 'red', main = 'LSA Position estimation error from real values', xlab = 'No. of iterations', ylab = 'log positioning error [m]')
lines(iter, log(abs(yy)), type = 'l', col = 'green')
lines(iter, log(abs(zz)), type = 'l', col = 'blue')

plot(iter, (abs(xx)), type = 'l', col = 'red', main = 'LSA Position estimation error from real values', xlab = 'No. of iterations', ylab = 'positioning error [m]')
lines(iter, (abs(yy)), type = 'l', col = 'green')
lines(iter, (abs(zz)), type = 'l', col = 'blue')


file.remove('razmakIteracija.txt','stvarnoOdstupanje.txt')