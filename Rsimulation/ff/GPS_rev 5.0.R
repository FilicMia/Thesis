rm(list=ls())

library(MASS)
library(matlib)
library(limSolve)
library(matrixcalc)

#options(error=recover)
npass <- 0 # po?etni postav brojila iteracija

S <- read.csv('satellites5.txt', header = FALSE) # u?itavanje polo?aja satelita

p <- read.csv('pseudoranges5a.txt', header = FALSE) # u?itavanje izmjerenih pseudoudaljenosti 
# datoteka xxxa.txt - nekorigirane pogre?ke
# datoteka xxxb.txt - korigirane pogre?ke

p <- as.vector(p[,1])
p <- t(t(p))

if(length(p) <= length(S)){
  stop('Not enough satellites - unable to estimate position. The script will quit.')
}

delta_x <- c(11,11,11,11) # postav po?etnih uvjeta za iteraciju - razlika susjednih iteracija (zaustavlja iteraciju)
x_0 <- t(t(c(0, 0, 0, 0))) # postav po?etnih uvjeta za iteraciju (procjena polo?aja)

# Definicija matrica i fizkalnih konstanti
J <- matrix(nrow = length(t(p)), ncol = length(S)+1)
R <- c(1,1)
dpr <- c(1,1)
W <- diag(length(t(p))) # matrica kovarijancija za Weighted LS solution - u po?etku postavljena kao jedini?na matrica
                         # za regular position estimation (pretpostavljena potpuna kompenzacija pogre?aka)
c <- 2.99792458E+08 # brzina svjetlosti [m/s], po GPS standardu
eps <- 1.0 # najve?a prihvatljiva pogre?ka komponente odre?ivanja polo?aja [m], eps = max(eps(x), eps(y), eps(z))

err <- c(0,0,0)
rlevel <- 11
start.time <- Sys.time()
while(eps < rlevel){ # iteracija - sve dok sve pogre?ke po komponentama ne budu manje od eps
for(i in 1:length(t(p))){
  R[i] <- sqrt((S[i,1] - x_0[1])^2 + (S[i,2] - x_0[2])^2 + (S[i,3] - x_0[3])^2) #udaljenost satelita i x_0
  
  dpr[i] <- p[i] - R[i] #b = pseudoudeljenost - izračunata udaljenost
  
  ssv <- 0
 
  for(j in 1:3){
    J[i,j] <- (x_0[j] - S[i,j])/R[i] #matrica parcijalnih derivacija
    ssv <- ssv + (S[i,j])^2
  }
  
  J[i, 4] <- c
  
  # Procjena kuta elevacije satelita
  d_x <- S[i,1] - x_0[1]
  d_y <- S[i,2] - x_0[2]
  d_z <- S[i,3] - x_0[3]
  
  
  A <- acos((S[i,1] * d_x + S[i,2] * d_y + S[i,3] * d_z)/ssv)
  ele <- pi/2 - A
  
  # [i, i]-ti element matrice kovarijancija 
  
  W[i, i] <- 1/(sin(ele))^2
  print(W)
  
  }

  # Postupak procjene polo?aja
  
  #dx <- qr.coef(qr(H), dpr) # Metoda A: non-WLSA
  
  dx <- chol2inv(t(J) %*% W %*% J) %*% t(J) %*% W %*% dpr # Metoda B: WLSA
  #D_T se uopće ne modelira...
  
  x_0 <- x_0 + dx
  
  npass <- npass + 1
  
  cat(c(npass, dx[1], dx[2], dx[3]),' \r',file="dx.txt", append=TRUE) # upisivanje vrijednosti dx radi kasnije analize brzine i to?nosti postupka
  err[1] <- x_0[1] - 918074.1038
  err[2] <- x_0[2] - 5703773.539
  err[3] <- x_0[3] -2693918.9285 
  cat(c(npass, err[1], err[2], err[3]),' \r',file="S.txt", append=TRUE)
  # Kontrola

  print(npass)
  print (S)
  print (dx)
  
  rlevel <- abs(dx[1]) + abs(dx[2]) + abs(dx[3])
  
}
closure <- Sys.time()

end.time <- Sys.time()

duration <- end.time - start.time

gr1 <- read.csv('dx.txt', header = FALSE, sep = '')
gr2 <- read.csv('S.txt', header = FALSE, sep = '')

np1 <- gr1$V1
d.x <- gr1$V2
d.y <- gr1$V3
d.z <- gr1$V4

plot(np1, log10(abs(d.x)), type = 'l', col = 'red', main = c('WLSA Time of execution in [s]=', round(duration, digits = 2)), xlab = 'No. of iterations', ylab = 'dx components')
lines(np1, log10(abs(d.y)), type = 'l', col = 'green')
lines(np1, log10(abs(d.z)), type = 'l', col = 'blue')

np1 <- gr2$V1
xx <- gr2$V2
yy <- gr2$V3
zz <- gr2$V4

plot(np1, log10(abs(xx)), type = 'l', col = 'red', main = 'WLSA Position estimation log10 error components', xlab = 'No. of iterations', ylab = 'positioning error [m]')
lines(np1, log10(abs(yy)), type = 'l', col = 'green')
lines(np1, log10(abs(zz)), type = 'l', col = 'blue')

file.remove('dx.txt','S.txt') # Brisanje privremenih daztoteka kori?tenih za crtanje dijagrama

### TO DO (kozmetika!)

# Izbor WLSA - non-WLSA obavlja se ru?no komentiranjem/ osloba?anjem L76 (non-WLSA) i L78 (WLSA)
# L112 i L121 ru?no popravljen naslov
# Ubaciti jo? jednu petlju u kojoj ?e prvo biti izveden jedan, a zatim drugi model (pomo?u signalne zastavice)