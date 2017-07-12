rm(list=ls())

library(MASS)
library(matlib)
library(limSolve)
library(matrixcalc)

iter = 0
niter = 100

c <- 2.99792458E+08 # brzina svjetlosti [m/s], po GPS standardu
R = read.csv('pseudoranges5a.txt', header = FALSE);
R <- as.matrix(p[,1])

#učitaj koordinate satelita
S = read.csv('satellites5.txt', header = FALSE)
S <- as.matrix(S)
nRows = dim(S)[1]
nCols = dim(S)[2]+1


if(length(p) < nCols){
  stop('Not enough satellites - unable to estimate position. The script will quit.')
}

delt <- c(11,11,11,11) # postav po?etnih uvjeta za iteraciju - razlika susjednih iteracija (zaustavlja iteraciju)
x_0 <- c(0, 0, 0, 0) # postav po?etnih uvjeta za iteraciju (procjena polo?aja)

# Definicija matrica i fizkalnih konstanti
J <- matrix(nrow = nRows, ncol = nCols)
b <- c(1,1)
W <- diag(nRows) # matrica kovarijancija za Weighted LS solution - u po?etku postavljena kao jedini?na matrica
# za regular position estimation (pretpostavljena potpuna kompenzacija pogre?aka)
c <- 2.99792458E+08 # brzina svjetlosti [m/s], po GPS standardu
eps <- 1.0 # najve?a prihvatljiva pogre?ka komponente odre?ivanja polo?aja [m], eps = max(eps(x), eps(y), eps(z))

rlevel <- 11
start.time <- Sys.time()
while(eps < rlevel){ # iteracija - sve dok sve pogre?ke po komponentama ne budu manje od eps
  for(i in 1:length(t(p))){
    # R[i] <- sqrt((S[i,1] - x_0[1])^2 + (S[i,2] - x_0[2])^2 + (S[i,3] - x_0[3])^2) #udaljenost satelita i x_0
    # 
    # dpr[i] <- p[i] - R[i] #b = pseudoudeljenost - izračunata udaljenost
    # 
    # ssv <- 0
    # 
    # for(j in 1:3){
    #   J[i,j] <- (x_0[j] - S[i,j])/R[i] #matrica parcijalnih derivacija
    #   ssv <- ssv + (S[i,j])^2
    # }
    
    #Matrica parcijalnihderivacija J of f, dijela
    
    J[i, 4] <- c
    
    # Procjena kuta elevacije satelita
    d_x <- S[i,1] - x_0[1]
    d_y <- S[i,2] - x_0[2]
    d_z <- S[i,3] - x_0[3]
    
    sr <- sqrt(ssv + (d_x^2 + d_y^2 + d_z^2))
    A <- acos((S[i,1] * d_x + S[i,2] * d_y + S[i,3] * d_z)/ssv)
    ele <- pi/2 - A
    
    # [i, i]-ti element matrice kovarijancija 
    
    W[i, i] <- 1/(sin(ele))^2
    
  }
  
  dx <- chol2inv(t(J) %*% W %*% J) %*% t(J) %*% W %*% dpr # Metoda B: WLSA
  #D_T se uopće ne modelira...
  
  x_0 <- x_0 + dx
  
  iter <- iter + 1
  
  cat(c(npass, dx[1], dx[2], dx[3]),' \r',file="razmakIteracija.txt", append=TRUE) # upisivanje vrijednosti dx radi kasnije analize brzine i to?nosti postupka
  err[1] <- x_0[1] - 918074.1038
  err[2] <- x_0[2] - 5703773.539
  err[3] <- x_0[3] -2693918.9285 
  cat(c(npass, err[1], err[2], err[3]),' \r',file="stvarnoOdstupanje.txt", append=TRUE)
  # Kontrola
  
  print(npass)
  print (S)
  print (dx)
  
  rlevel <- abs(dx[1]) + abs(dx[2]) + abs(dx[3])
  
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