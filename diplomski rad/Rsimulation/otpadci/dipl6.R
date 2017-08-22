rm(list=ls())

library(MASS)
library(matlib)
library(limSolve)
library(matrixcalc) #NISTA

iter = 0
niter = 100
#option = 1 # W se bira kao dijaginalna matrica s vrijednostima 1/sin(Ei) na dijaginali
#option = 2 # W se bira kao dijaginalna matrica s vrijednostima a_ele**2+b_ele**2/sin(Ei) na dijaginali
#option = 3 # W se bira kao dijaginalna matrica s vrijednostima a**2/sin(Ei+psi) na dijaginali
option = 3 #<- ok, izbjegava se singularitet u 0
solution = "Ch" #ili "Ch"; odabir modela rješavanja sustava.
doChange = TRUE

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
x_0 <- c(11, 11, 11, 0) #nije (0,0,0,0) kako bi se izbjegnuo singularitet prilikom prve iteracije i izračunavanja
# kuta elevacije

unutar = append(S,rep(-c,nRows))
SC = matrix(unutar,nRows,nCols) # [x_i,y_i,z_i,c] 

start.time <- Sys.time()# za regular position estimation (pretpostavljena potpuna kompenzacija pogre?aka)
eps <- 1 #eps > max(delt(x), delt(y), delt(z)), kriterij zaustavljanja
#zanimljive stvari se dogadaju za 10^3.2

A_iter <- matrix(nrow = nRows, ncol=nCols)

rlevel <- 10^5 #max(delt(x), delt(y), delt(z))
b <- R

#ako pomak pojedine iteracije padne ispod određene vrijednosti, ne mijenjamo više tu koordinatu
change <- rep(TRUE,nCols-1) #indikatori mijenjanja koordinata (x,y,z) redom.
chl = FALSE
start.time <- Sys.time() 

while(eps < rlevel){ 
  # iteracija - sve dok sve pogreske po komponentama ne budu manje od eps
  x_iter = c(x_0[1:3],0) # samo c , a ne c-d_T
  AA = t(apply(SC, 1, function(x) (x_iter - x))) #zbog lakse derivacije je x_-x
  D = sqrt((AA*AA)%*%c(1,1,1,0))
  DD = matrix(append(rep(D,3),rep(1,nRows)),nRows,nCols)
  
  A_iter = AA/DD #J_k
  
  #Procjena kuta elevacije satelita
  #n = (x_iter[i,1],x_iter[i,2],x_iter[i,3])
  #s = (S[i,1],S[i,2], S[i,3] )
  xyz.coords.matrix = matrix(rep(x_iter[1:nCols-1],nRows),nRows,nCols-1)
  
  xyz.coords.ssv <- sqrt(xyz.coords.matrix**2%*%c(1,1,1))#zbroj svih kooordinata na kvadrat
  xyz.coords.ssv = matrix(rep(xyz.coords.ssv,3),nRows,nCols-1)
  
  function.ele <- function(angle){
    if(angle < 0){
      a = -angle + pi/2
      return(pi-a)
    }else{
      return(pi/2 - angle)
    }
  }
  
  if(option==1){
    #budući da uvijek računamo x_iter - x, treba se nanovo izračunati x - x_iter. Smjer
    #vektora je obrnut. -AA
    
    D_xyz = -AA[1:nRows,1:(nCols-1)] # zraka x_iter do Si
    D_xyz.ssv <- sqrt(D_xyz**2%*%c(1,1,1))
    D_xyz.ssv = matrix(rep(D_xyz.ssv,3),nRows,nCols-1)
    
    E <- acos(((xyz.coords.matrix/xyz.coords.ssv)*(D_xyz/D_xyz.ssv))%*%c(1,1,1))#kut izmedu normale tangencijalne ravnine i zrake prijemnik-satelit.
    ele <- apply(E,1,function.ele)
    
    Wii = 1/(tan(ele))^2 
    print("ele")
    print(ele)
    
  }else if(option==2){
    D_xyz = -A_iter[1:nRows,1:(nCols-1)]
    D_xyz.ssv <- sqrt(D_xyz**2%*%c(1,1,1))
    D_xyz.ssv = matrix(rep(D_xyz.ssv,3),nRows,nCols-1)
    
    E <- acos(((xyz.coords.matrix/xyz.coords.ssv)*(D_xyz/D_xyz.ssv))%*%c(1,1,1))
    ele <- apply(E,1,function.ele)
    
    Wii = a_ele**2+b_ele**2/(tan(ele))^2
    print("ele")
    print(ele)
  }else if(option == 3){ 
    D_xyz = -A_iter[1:nRows,1:(nCols-1)]
    D_xyz.ssv <- sqrt(D_xyz**2%*%c(1,1,1))
    D_xyz.ssv = matrix(rep(D_xyz.ssv,3),nRows,nCols-1)
    
    E <- acos(((xyz.coords.matrix/xyz.coords.ssv)*(D_xyz/D_xyz.ssv))%*%c(1,1,1))
    ele <- apply(E,1,function.ele)
    
    Wii = a**2/((tan(ele+0.01)+psi)^2)
    print("ele")
    print(ele)
    
  }
  
  W = diag(Wii)
  b <- R[,1] - D - c*x_0[nCols] #x_0[4] modelira d_T i delta = delta_{d_T}
  
  if(solution == "Ch"){
    dx <- chol2inv(t(A_iter) %*% W %*% A_iter) %*% t(A_iter) %*% W %*% b
  } else if (solution == "QR") {
    #drugikorijen iz W
    dx <- qr.coef(qr(sqrt(W)%*%A_iter), sqrt(W)%*%b)
  } else if (solution == "ChQR"){
    if(10^2 > rlevel){ #Ch do određene točnosti, inicijalna vrijednost
      dx <- qr.coef(qr(sqrt(W)%*%A_iter), sqrt(W)%*%b)
      solution = "Ch"
    } else {
      #Za option = 3 i eps = 10^3.2 dobivamo da pogreška x koordinate pada ispod 10^2 što je odlično.
      dx <- chol2inv(t(A_iter) %*% W %*% A_iter) %*% t(A_iter) %*% W %*% b
    } 
  } else {
    # } else if (solution == "ChQR10_2"){
    #   if(dx[1] > 10^3){ #Ch do određene točnosti
    #     dx <- chol2inv(t(A_iter) %*% W %*% A_iter) %*% t(A_iter) %*% W %*% b
    #   } else {
    #     dx <- qr.coef(qr(sqrt(W)%*%A_iter), sqrt(W)%*%b)
    #     dx <- svd.inverse(t(A_iter) %*% W %*% A_iter) %*% t(A_iter) %*% W %*% b
    #     cat(c(iter, err[1], err[2], err[3]),' \r',file="manjiDelta.txt", append=TRUE)
    #     solution <- "QR"
    #   }
    # } else {
    print("Nepoznata opcija rješavanja sustava jednadžbi.")
    exit(-1)
  }
  #dx <- svd.inverse(t(A_iter) %*% W %*% A_iter) %*% t(A_iter) %*% W %*% b # najpreciznija
  
  x_0 <- x_0 + dx*c(change,TRUE)
  
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
  print(change)
  
  if(doChange){
    change <- abs(dx[1:3,1]) > eps #Je li iteracija zaustavlja po komponentama. Kada određena komp od delta_x
    # postane dovoljno mala, njezina komponenta od x_iter se više ne mijenja.
  }
  #dovodi do poboljšanja ukoliko je broj iteracija veći
  rlevel <- max(abs(dx[1:3,1]))
  
}


end.time <- Sys.time()

timediff <- end.time - start.time

d_iter <- read.csv('razmakIteracija.txt', header = FALSE, sep = '')
err <- read.csv('stvarnoOdstupanje.txt', header = FALSE, sep = '')

iter <- d_iter$V1
d.x <- d_iter$V2
d.y <- d_iter$V3
d.z <- d_iter$V4

plot(iter, log10(abs(d.x)), type = 'l', col = 'red', main = c('LSA vrijeme izvršavanja u  [s]=', round(timediff, digits = 2)), xlab = 'broj iteracija', ylab = 'log10 deltaX komponenti')
lines(iter, log10(abs(d.y)), type = 'l', col = 'green')
lines(iter, log10(abs(d.z)), type = 'l', col = 'blue')
legend(
  "topright", legend=c("x","y", "z"), col=c("red", "green", "blue"),   pch=15
)

plot(iter, (abs(d.x)), type = 'l', col = 'red', main = c('LSA vrijeme izvršavanja u  [s]=', round(timediff, digits = 2)), xlab = 'broj iteracija', ylab = 'deltaX komponente')
lines(iter, (abs(d.y)), type = 'l', col = 'green')
lines(iter, (abs(d.z)), type = 'l', col = 'blue')
legend(
  "topright", legend=c("x","y", "z"), col=c("red", "green", "blue"),   pch=15
)

iter <- err$V1
xx <- err$V2
yy <- err$V3
zz <- err$V4

plot(iter, ylim=c(1.5,5),log10(abs(xx)), type = 'l', col = 'red', main = 'LSA odstupanje procjene položaja \n od stvarne vrijednosti', xlab = 'broj iteracija', ylab = 'log10 odstupanje [m]')
lines(iter, log10(abs(yy)), type = 'l', col = 'green')
lines(iter, log10(abs(zz)), type = 'l', col = 'blue')
legend(
  "topleft", legend=c("x","y", "z"), col=c("red", "green", "blue"),   pch=15
)

plot(iter, (abs(xx)), type = 'l', col = 'red', main = 'LSA odstupanje procjene položaja \n od stvarne vrijednosti', xlab = 'broj iteracija', ylab = 'odstupanje [m]')
lines(iter, (abs(yy)), type = 'l', col = 'green')
lines(iter, (abs(zz)), type = 'l', col = 'blue')
legend(
  "topleft", legend=c("x","y", "z"), col=c("red", "green", "blue"),   pch=15
)

file.remove('razmakIteracija.txt','stvarnoOdstupanje.txt')
#
#Sustav izgleda nadasve nestabilan, i ukoliko se koristibolja metoda računanja inverza 
# i pronalaska rješenja,
# algoritam konvergira, ali k krivom rješenju.
# Sporija konvergencija uzrokovana metodom Choleskoga, pomaže da algoritam konvergita ka 
# rješenju veće točnosti.
# Zaključujemo kako je model potrebno doraditi kako bi se postigla 
# 
# rješenju.