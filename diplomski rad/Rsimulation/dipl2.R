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
delt = c(1000,1000,1000,3)
nRows = dim(S)[1]
nCols = dim(S)[2]+1

realPosition = c(918074.1038,5703773.539,2693918.9285,0)

inside = append(S,rep(-c,nRows))
RS = matrix(inside,nRows,nCols) # [x_i,y_i,z_i,d_i] 

iter = 0
niter = 1000

err <- c(11,11,11,11)
#while(norm(t(delt)) > 1){
start.time <- Sys.time()

b = R
while(iter < 100*niter && (max(abs(delt[1:3])) > 1000 || iter == 0)){

  x_iter = c(x_0[1:3],0) #delta x_iter
  AA = t(apply(RS, 1, function(x) (x_iter - x))) 
  
  #udaljenost satelita od procjenjenog položaja
  D = sqrt(AA**2%*%c(1,1,1,0))
  DD = matrix(append(rep(D,3),rep(1,nRows)),nRows,nCols) 
  
  A_iter = AA/DD
  #https://www.math.ucla.edu/~anderson/rw1001/library/base/html/qr.html
  #exaple -> zadnji red
  #QR = qr(A, tol = 1e-07 , LAPACK = TRUE)
  #R = qr.R(QR,complete=TRUE)
  #Q = qr.Q(QR,complete=TRUE) #vraca N*N matricu <- cijelu QR matricu
  
  # rješava sustav Ax=b koristeći QR faktorizaciju
  delt <- qr.coef(qr(A_iter), b) 
  
  x_0 = x_0 + delt #(x,y,z,dT)
  b = R - D - c*x_0[nCols]
  
  #upisivanje vrijednosti dx radi kasnije analize brzine i točnosti postupka
  if(iter%%10 == 0){
    cat(c(iter, delt[1:3]),' \r',file="razmakIteracija.txt", append=TRUE) 
    err <- x_0 - realPosition
    cat(c(iter, err[1:3]),' \r',file="stvarnoOdstupanje.txt", append=TRUE) 
    print(A_iter) 
  }
  if(abs(delt[1]) < 100 ||
     abs(delt[2]) < 100 ||
     abs(delt[3]) < 100){
    cat(c(iter, delt[1:3]),' \r',file="razmakIteracija100.txt", append=TRUE) 
    cat(c(iter, err[1:3]),' \r',file="stvarnoOdstupanje100.txt", append=TRUE) 
  }
  
  if(abs(delt[1]) < 1000 ||
     abs(delt[2]) < 1000 ||
     abs(delt[3]) < 1000 ){
    cat(c(iter, delt[1:3]),' \r',file="razmakIteracija1000.txt", append=TRUE) 
    cat(c(iter, err[1:3]),' \r',file="stvarnoOdstupanje1000.txt", append=TRUE) 
  }
  
  iter = iter +1
 
}

end.time <- Sys.time()
timediff <- end.time - start.time

#grafički prikaz preciznosti i točnosti odredivanja položaja kroz iteracije
d_iter <- read.csv('razmakIteracija.txt', header = FALSE, sep = '')
err <- read.csv('stvarnoOdstupanje.txt', header = FALSE, sep = '')

iter <- d_iter$V1
d.x <- d_iter$V2
d.y <- d_iter$V3
d.z <- d_iter$V4

plot(iter, log10(abs(d.x)), type = 'l', col = 'red', main = c('LSA vrijeme izvršavanja u [s]=', round(timediff, digits = 2)), xlab = 'Broj iteracija', ylab = 'log10 delta_x po komponentama')
lines(iter, log10(abs(d.y)), type = 'l', col = 'green')
lines(iter, log10(abs(d.z)), type = 'l', col = 'blue')

plot(iter, (abs(d.x)), type = 'l', col = 'red', main = c('LSA vrijeme izvršavanja u [s]=', round(timediff, digits = 2)), xlab = 'Broj iteracija', ylab = 'delta_x po komponentama')
lines(iter, (abs(d.y)), type = 'l', col = 'green')
lines(iter, (abs(d.z)), type = 'l', col = 'blue')

iter <- err$V1
xx <- err$V2
yy <- err$V3
zz <- err$V4

plot(iter, log10(abs(xx)), type = 'l', col = 'red', main = 'LSA pogreška od stvarnog položaja', xlab = 'Broj iteracija', ylab = 'log10 pogreška [m]')
lines(iter, log10(abs(yy)), type = 'l', col = 'green')
lines(iter, log10(abs(zz)), type = 'l', col = 'blue')

plot(iter, (abs(xx)), type = 'l', col = 'red', main = 'LSA pogreška od stvarnog položaja', xlab = 'Broj iteracija', ylab = 'pogreška [m]')
lines(iter, (abs(yy)), type = 'l', col = 'green')
lines(iter, (abs(zz)), type = 'l', col = 'blue')


file.remove('razmakIteracija.txt','stvarnoOdstupanje.txt')

d_iter100 <- read.csv('razmakIteracija100.txt', header = FALSE, sep = '')
err100 <- read.csv('stvarnoOdstupanje100.txt', header = FALSE, sep = '')

iter <- d_iter100$V1
d.x <- d_iter100$V2
d.y <- d_iter100$V3
d.z <- d_iter100$V4

plot(iter, log10(abs(d.x)), type = 'p', col = 'red', 
     main = c('LSA vrijeme izvršavanja u [s]=', round(timediff, digits = 2)), xlab = 'Broj iteracija', ylab = 'log10 delta_x po komponentama')
lines(iter, log10(abs(d.y)), type = 'p', col = 'green')
lines(iter, log10(abs(d.z)), type = 'p', col = 'blue')

plot(iter, (abs(d.x)), type = 'p', col = 'red', 
     main = c('(< 100 [m] bar 1 koordinati) LSA vrijeme izvršavanja u [s]=', round(timediff, digits = 2)), xlab = 'Broj iteracija', ylab = 'delta_x po komponentama')
lines(iter, (abs(d.y)), type = 'p', col = 'green')
lines(iter, (abs(d.z)), type = 'p', col = 'blue')

iter <- err100$V1
xx <- err100$V2
yy <- err100$V3
zz <- err100$V4

plot(iter, log10(abs(xx)), type = 'p', col = 'red', 
     main = 'LSA pogreška od stvarnog položaja \n (< 100 [m] bar 1 koordinati)', xlab = 'Broj iteracija', ylab = 'log10 pogreška [m]')
lines(iter, log10(abs(yy)), type = 'p', col = 'green')
lines(iter, log10(abs(zz)), type = 'p', col = 'blue')

plot(iter, (abs(xx)), type = 'p', col = 'red', 
     main = 'LSA pogreška od stvarnog položaja \n (< 100 [m] bar 1 koordinati)', xlab = 'Broj iteracija', ylab = 'pogreška [m]')
lines(iter, (abs(yy)), type = 'p', col = 'green')
lines(iter, (abs(zz)), type = 'p', col = 'blue')

# < 1000 [m] barem 2 koordinate
d_iter1000 <- read.csv('razmakIteracija1000.txt', header = FALSE, sep = '')
err1000 <- read.csv('stvarnoOdstupanje1000.txt', header = FALSE, sep = '')

iter <- d_iter1000$V1
d.x <- d_iter1000$V2
d.y <- d_iter1000$V3
d.z <- d_iter1000$V4

plot(iter, log10(abs(d.x)), type = 'p', col = 'red', 
     main = c('(1000) LSA vrijeme izvršavanja u [s]=', round(timediff, digits = 2)), xlab = 'Broj iteracija', ylab = 'log10 delta_x po komponentama')
lines(iter, log10(abs(d.y)), type = 'p', col = 'green')
lines(iter, log10(abs(d.z)), type = 'p', col = 'blue')

plot(iter, (abs(d.x)), type = 'p', col = 'red', 
     main = c('(< 1000 [m] bar 1 koordinati)\n LSA vrijeme izvršavanja u [s]=', round(timediff, digits = 2)), xlab = 'Broj iteracija', ylab = 'delta_x po komponentama')
lines(iter, (abs(d.y)), type = 'p', col = 'green')
lines(iter, (abs(d.z)), type = 'p', col = 'blue')

iter <- err1000$V1
xx <- err1000$V2
yy <- err1000$V3
zz <- err1000$V4

plot(iter, log10(abs(xx)), type = 'p', col = 'red', 
     main = 'LSA pogreška od stvarnog položaja \n (< 1000 [m] bar 1 koordinati)', xlab = 'Broj iteracija', ylab = 'log10 pogreška [m]')
lines(iter, log10(abs(yy)), type = 'p', col = 'green')
lines(iter, log10(abs(zz)), type = 'p', col = 'blue')

plot(iter, (abs(xx)), type = 'p', col = 'red', 
     main = 'LSA pogreška od stvarnog položaja \n (< 1000 [m] bar 1 koordinati)', xlab = 'Broj iteracija', ylab = 'pogreška [m]')
lines(iter, (abs(yy)), type = 'p', col = 'green')
lines(iter, (abs(zz)), type = 'p', col = 'blue')

file.remove('razmakIteracija100.txt','stvarnoOdstupanje100.txt','razmakIteracija1000.txt','stvarnoOdstupanje1000.txt')