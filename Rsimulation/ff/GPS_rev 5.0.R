rm(list=ls())

start.time <- Sys.time()

library(MASS)
library(matlib)
library(limSolve)
library(matrixcalc)

#options(error=recover)
getwd()
npass <- 0 # po?etni postav brojila iteracija

x <- read.csv('satellites5.txt', header = FALSE) # u?itavanje polo?aja satelita

pr <- read.csv('pseudoranges5a.txt', header = FALSE) # u?itavanje izmjerenih pseudoudaljenosti 
# datoteka xxxa.txt - nekorigirane pogre?ke
# datoteka xxxb.txt - korigirane pogre?ke

pr <- as.vector(pr[,1])
pr <- t(t(pr))

if(length(pr) < length(x)){
  stop('Not enough satellites - unable to estimate position. The script will quit.')
}

dx <- c(11,11,11,11) # postav po?etnih uvjeta za iteraciju - razlika susjednih iteracija (zaustavlja iteraciju)
X <- t(t(c(0, 0, 0, 0))) # postav po?etnih uvjeta za iteraciju (procjena polo?aja)

# Definicija matrica i fizkalnih konstanti
H <- matrix(nrow = length(t(pr)), ncol = length(x)+1)
R <- c(1,1)
dpr <- c(1,1)
C <- diag(length(t(pr))) # matrica kovarijancija za Weighted LS solution - u po?etku postavljena kao jedini?na matrica
                         # za regular position estimation (pretpostavljena potpuna kompenzacija pogre?aka)
c <- 2.99792458E+08 # brzina svjetlosti [m/s], po GPS standardu
eps <- 1.0 # najve?a prihvatljiva pogre?ka komponente odre?ivanja polo?aja [m], eps = max(eps(x), eps(y), eps(z))

err <- c(0,0,0)
rlevel <- 11

while(eps < rlevel){ # iteracija - sve dok sve pogre?ke po komponentama ne budu manje od eps
for(i in 1:length(t(pr))){
  R[i] <- sqrt((x[i,1] - X[1])^2 + (x[i,2] - X[2])^2 + (x[i,3] - X[3])^2)
  
  dpr[i] <- pr[i] - R[i]
  
  ssv <- 0
 
  for(j in 1:3){
    H[i,j] <- (X[j] - x[i,j])/R[i]
    ssv <- ssv + (x[i,j])^2
  }
  
  H[i, 4] <- c
  
  # Procjena kuta elevacije satelita
  d_x <- x[i,1] - X[1]
  d_y <- x[i,2] - X[2]
  d_z <- x[i,3] - X[3]
  
  sr <- sqrt(ssv + (d_x^2 + d_y^2 + d_z^2))
  A <- acos((x[i,1] * d_x + x[i,2] * d_y + x[i,3] * d_z)/ssv)
  ele <- pi/2 - A
  
  # [i, i]-ti element matrice kovarijancija 
  
  C[i, i] <- 1/(sin(ele))^2
  
  }

  # Postupak procjene polo?aja
  
  #dx <- qr.coef(qr(H), dpr) # Metoda A: non-WLSA
  
  dx <- chol2inv(t(H) %*% C %*% H) %*% t(H) %*% C %*% dpr # Metoda B: WLSA
  
  X <- X + dx
  
  npass <- npass + 1
  
  cat(c(npass, dx[1], dx[2], dx[3]),' \r',file="dx.txt", append=TRUE) # upisivanje vrijednosti dx radi kasnije analize brzine i to?nosti postupka
  err[1] <- X[1] - 918074.1038
  err[2] <- X[2] - 5703773.539
  err[3] <- X[3] -2693918.9285 
  cat(c(npass, err[1], err[2], err[3]),' \r',file="X.txt", append=TRUE)
  # Kontrola

  print(npass)
  print (X)
  print (dx)
  
  rlevel <- abs(dx[1]) + abs(dx[2]) + abs(dx[3])
  
}
closure <- Sys.time()

end.time <- Sys.time()

duration <- end.time - start.time

gr1 <- read.csv('dx.txt', header = FALSE, sep = '')
gr2 <- read.csv('X.txt', header = FALSE, sep = '')

np1 <- gr1$V1
d.x <- gr1$V2
d.y <- gr1$V3
d.z <- gr1$V4

plot(np1, log(abs(d.x)), type = 'l', col = 'red', main = c('WLSA Time of execution in [s]=', round(duration, digits = 2)), xlab = 'No. of iterations', ylab = 'dx components')
lines(np1, log(abs(d.y)), type = 'l', col = 'green')
lines(np1, log(abs(d.z)), type = 'l', col = 'blue')

np1 <- gr2$V1
xx <- gr2$V2
yy <- gr2$V3
zz <- gr2$V4

plot(np1, log(abs(xx)), type = 'l', col = 'red', main = 'WLSA Position estimation error components', xlab = 'No. of iterations', ylab = 'positioning error [m]')
lines(np1, log(abs(yy)), type = 'l', col = 'green')
lines(np1, log(abs(zz)), type = 'l', col = 'blue')

file.remove('dx.txt','X.txt') # Brisanje privremenih daztoteka kori?tenih za crtanje dijagrama

### TO DO (kozmetika!)

# Izbor WLSA - non-WLSA obavlja se ru?no komentiranjem/ osloba?anjem L76 (non-WLSA) i L78 (WLSA)
# L112 i L121 ru?no popravljen naslov
# Ubaciti jo? jednu petlju u kojoj ?e prvo biti izveden jedan, a zatim drugi model (pomo?u signalne zastavice)