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
nRows = dim(S)[1]
nCols = dim(S)[2]+1

realPosition = c(918074.1038,5703773.539,2693918.9285,0)

unutar = append(S,rep(-c,nRows))
RS = matrix(unutar,nRows,nCols) # [x_i,y_i,z_i,d_i] 

iter = 0
niter = 1000

err <- c(11,11,11,11)
#while(norm(t(delt)) > 1){
start.time <- Sys.time()

b = R
while(iter < niter){
  
  x_iter = c(x_0[1:3],0) # samo c , a ne c-d_T
  AA = t(apply(RS, 1, function(x) (x_iter - x))) #zbog lakše derivacije je x_-x
  D = sqrt(AA**2%*%c(1,1,1,0))
  DD = matrix(append(rep(D,3),rep(1,nRows)),nRows,nCols)
  
  
  A_iter = AA/DD
  delt <- qr.coef(qr(A_iter), b) # rješava sustav Ax=b koristeći QR faktorizaciju
  
  x_0 = x_0 + delt #(x,y,z,c*dT)
  
  b = R - D
  if(iter%%10 == 0){
    cat(c(iter, delt[1:3]),' \r',file="razmakIteracija.txt", append=TRUE) # upisivanje vrijednosti dx radi kasnije analize brzine i to?nosti postupka
    err <- x_0 - realPosition
    cat(c(iter, err[1:3]),' \r',file="stvarnoOdstupanje.txt", append=TRUE) 
    print(A_iter) 
  }
  
  iter = iter +1
  
}

end.time <- Sys.time()

timediff <- end.time - start.time

d_iter <- read.csv('razmakIteracija.txt', header = FALSE, sep = '')
err <- read.csv('stvarnoOdstupanje.txt', header = FALSE, sep = '')

iter <- d_iter$V1
d.x <- d_iter$V2
d.y <- d_iter$V3
d.z <- d_iter$V4

plot(iter, log10(abs(d.x)), type = 'l', col = 'red', main = c('LSA Time of execution in [s]=', round(timediff, digits = 2)), xlab = 'No. of iterations', ylab = 'log10 of d_X components')
lines(iter, log10(abs(d.y)), type = 'l', col = 'green')
lines(iter, log10(abs(d.z)), type = 'l', col = 'blue')

plot(iter, (abs(d.x)), type = 'l', col = 'red', main = c('LSA Time of execution in [s]=', round(timediff, digits = 2)), xlab = 'No. of iterations', ylab = 'd_X components')
lines(iter, (abs(d.y)), type = 'l', col = 'green')
lines(iter, (abs(d.z)), type = 'l', col = 'blue')

iter <- err$V1
xx <- err$V2
yy <- err$V3
zz <- err$V4

plot(iter, log10(abs(xx)), type = 'l', col = 'red', main = 'LSA Position estimation error from real values', xlab = 'No. of iterations', ylab = 'log10 positioning error [m]')
lines(iter, log10(abs(yy)), type = 'l', col = 'green')
lines(iter, log10(abs(zz)), type = 'l', col = 'blue')

plot(iter, (abs(xx)), type = 'l', col = 'red', main = 'LSA Position estimation error from real values', xlab = 'No. of iterations', ylab = 'positioning error [m]')
lines(iter, (abs(yy)), type = 'l', col = 'green')
lines(iter, (abs(zz)), type = 'l', col = 'blue')


file.remove('razmakIteracija.txt','stvarnoOdstupanje.txt')