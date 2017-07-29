rm(list=ls())
#ubrzanje, 0.2 naprema 0.43

library(MASS)
library(matlib)
library(limSolve)
library(matrixcalc)

iter = 0
niter = 100

c <- 2.99792458E+08 # brzina svjetlosti [m/s], po GPS standardu
R = read.csv('pseudoranges5a.txt', header = FALSE);
R <- as.matrix(R[,1])
realPosition <- c(918074.1038,5703773.539,2693918.9285)

#ucitaj koordinate satelita
S = read.csv('satellites5.txt', header = FALSE)
S <- as.matrix(S)
nRows = dim(S)[1]
nCols = dim(S)[2]+1

#konstante iteracije
RR <- rep(0,nRows)
dpr <- c(1,1)
W <- diag(nRows)
#option = 1 # W se bira kao dijaginalna matrica s vrijednostima 1/sin(Ei) na dijaginali
#option = 2 # W se bira kao dijaginalna matrica s vrijednostima a_ele**2+b_ele**2/sin(Ei) na dijaginali
#option = 3 # W se bira kao dijaginalna matrica s vrijednostima a**2/sin(Ei+psi) na dijaginali
option = 3 #<- ok, izbjegava se singularitet u 0
a_ele = 1
b_ele= 2

a = 1
psi = 0.5


if(dim(R)[1] < nCols){
  stop('Not enough satellites - unable to estimate position. The script will quit.')
}

#pocetni uvjeti
delt <- c(11,11,11,11)
err <- c(11,11,11,11)
x_0 <- c(0, 0, 0, 0) 

unutar = append(S,rep(-c,nRows))
SC = matrix(unutar,nRows,nCols) # [x_i,y_i,z_i,c] 

start.time <- Sys.time()# za regular position estimation (pretpostavljena potpuna kompenzacija pogre?aka)
eps <- 1.0 #eps > max(delt(x), delt(y), delt(z)), kriterij zaustavljanja
A_iter <- matrix(nrow = nRows, ncol=nCols)

rlevel <- 11 #max(delt(x), delt(y), delt(z))
b <- R
start.time <- Sys.time() #while(norm(t(delt)) > 1){
while(eps < rlevel ){ 
  # iteracija - sve dok sve pogreske po komponentama ne budu manje od eps
  x_iter = c(x_0[1:3],0) # samo c , a ne c-d_T
  AA = t(apply(SC, 1, function(x) (x_iter - x))) #zbog lakse derivacije je x_-x
  D = sqrt((AA*AA)%*%c(1,1,1,0))
  DD = matrix(append(rep(D,3),rep(1,nRows)),nRows,nCols)
    
  A_iter = AA/DD #J_k

  #Procjena kuta elevacije satelita
  D_xyz = AA[1:nRows,1:(nCols-1)] # zraka x_iter do Si
  #n = (S[i,1],S[i,2],S[i,3])
  #s = (S[i,1],S[i,2], 0 )
  
  ssv <- S**2%*%c(1,1,1)#zbroj svih kooordinata satelita na kvadrat
  ssv = matrix(rep(ssv,3),nRows,nCols-1)
    
  E <- acos(((S*D_xyz)/ssv)%*%c(1,1,1)) #?? kaj nije da se dijeli i s nornom od D_xyz??
  ele <- pi/2 - E
  
  if(option==2){
    D_xyz = A_iter[1:nRows,1:(nCols-1)]
    E <- acos(((S*D_xyz)/ssv)%*%c(1,1,1)) #?? kaj nije da se dijeli i s nornom od D_xyz??
    ele <- pi/2 - E
    
    Wii = a_ele**2+b_ele**2/(sin(ele))^2
    print("ele")
    print(ele)
  }else if(option == 3){ 
    D_xyz = A_iter[1:nRows,1:(nCols-1)]
    E <- acos(((S*D_xyz)/ssv)%*%c(1,1,1)) #?? kaj nije da se dijeli i s nornom od D_xyz??
    ele <- pi/2 - E
    
    Wii = a**2/(sin(ele+psi))^2
    print("ele")
    print(ele)
    
  }else{
    Wii = 1/(sin(ele))^2 
    print("ele")
    print(ele)
  }
  W = diag(Wii[,1])
  b <- R[,1] - D - c*x_0[nCols] #x_0[4] modelira d_T i delta = delta_{d_T}
    
  dx <- chol2inv(t(A_iter) %*% W %*% A_iter) %*% t(A_iter) %*% W %*% b # 
  #dx <- svd.inverse(t(A_iter) %*% W %*% A_iter) %*% t(A_iter) %*% W %*% b # najpreciznija
  
  #drugikorijen iz W
  #dx <- qr.coef(qr(sqrt(W)%*%A_iter), sqrt(W)%*%b)
  
  x_0 <- x_0 + dx
  
  iter <- iter + 1
  
  cat(c(iter, dx[1], dx[2], dx[3]),' \r',file="razmakIteracija.txt", append=TRUE) # upisivanje vrijednosti dx radi kasnije analize brzine i to?nosti postupka
  err[1] <- x_0[1] - 918074.1038
  err[2] <- x_0[2] - 5703773.539
  err[3] <- x_0[3] -2693918.9285 
  cat(c(iter, err[1], err[2], err[3]),' \r',file="stvarnoOdstupanje.txt", append=TRUE)
  # Kontrola
  
  print(iter)
  print (S)
  print (dx)
  
  rlevel <- max(abs(dx[1]), abs(dx[2]), abs(dx[3]))
  
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
#
#Sustav izgleda nadasve nestabilan, i ukoliko se koristibolja metoda računanja inverza 
# i pronalaska rješenja,
# algoritam konvergira, ali k krivom rješenju.
# Sporija konvergencija uzrokovana metodom Choleskoga, pomaže da algoritam konvergita ka 
# rješenju veće točnosti.
# Zaključujemo kako je model potrebno doratiti kako bi se postigla 
# 
# rješenju.

#dipl1 je najbolja. !!!