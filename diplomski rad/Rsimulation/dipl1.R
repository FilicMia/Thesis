library("MASS","matrixcalc")
#pseudo-udaljenosti
c <- 2.99792458E+08 # brzina svjetlosti [m/s], po GPS standardu
R = read.csv('pseudorangesb.txt', header = FALSE);
R <- as.matrix(R[,1])
#učitaj koordinate satelita
S = read.csv('satellites.txt', header = FALSE)
S <- as.matrix(S)
nRows = dim(S)[1]
nCols = dim(S)[2]+1

#konstante iteracije
W <- diag(nRows)

a_ele = 1
b_ele= 2

a = 1
psi = 0.5

change <- rep(TRUE,nCols-1) #indikatori mijenjanja koordinata (x,y,z) redom.
x_0 = c(1,1,1,1) #[x,y,z,d_T] d_t se kasnije množi sa c da bi se oduzeo od [x_i,y_i,z_i,d] 
delt = c(3,3,3,3)

tezine = FALSE#ako želimo algoritam težinskih najmanjih kvadrata
#option = 1 # W se bira kao dijaginalna matrica s vrijednostima 1/sin(Ei) na dijaginali
#option = 2 # W se bira kao dijaginalna matrica s vrijednostima a_ele**2+b_ele**2/sin(Ei) na dijaginali
#option = 3 # W se bira kao dijaginalna matrica s vrijednostima a**2/sin(Ei+psi) na dijaginali
option = 3 #<- ok, izbjegava se singularitet u 0
solution = "Ch" #ili "Ch"; odabir modela rješavanja sustava.
elevation = 1 #način računanja kuta elevacije je jedna 1 za satelit-prijemnik i satelit-ishodište
#2 za prijemnik-satelit i normala tang. ravnine u točki prijemnika na sferu radijusa udaljenosti prijemnika od ishodišta
#slike rezultata spremnjene s nastavkom new.

realPosition = c(918074.1038,5703773.539,2693918.9285,0)

unutar = append(S,R)
RS = matrix(unutar,dim(S)[1],dim(S)[2]+1) # [x_i,y_i,z_i,d] 

iter = 0
niter = 100
err <- c(11,11,11,11)
#while(norm(t(delt)) > 1){
start.time <- Sys.time()

while(iter < 100*niter && (max(abs(delt[1:3])) > 1 || iter == 0)){
  x_ = c(x_0[1:3],x_0[4]*c)
  P = t(apply(RS, 1, function(x) (x-x_))) #x-x_0
  D = sqrt((P*P)%*%c(1,1,1,0))
  DD = matrix(append(rep(D,3),rep(1,nRows)),nRows,nCols)
  
  PI = P
  PI[,dim(S)[2]+1] = - c*P[,dim(S)[2]+1] #4. red s -c
  
  PX = P*P
  b = PX%*%c(-1,-1,-1,1)
  A_iter = -2*PI
  
  if(tezine){
    #Procjena kuta elevacije satelita
    #n = (x_iter[i,1],x_iter[i,2],x_iter[i,3])
    #s = (S[i,1],S[i,2], S[i,3] )
    if(elevation == 1){
      xyz.coords.matrix = matrix(S,nRows,nCols-1)
    } else {
      xyz.coords.matrix = matrix(rep(x_iter[1:nCols-1],nRows),nRows,nCols-1)
      print(2)
    }
    
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
      
      D_xyz = -P[1:nRows,1:(nCols-1)] # zraka x_iter do Si
      D_xyz.ssv <- sqrt(D_xyz**2%*%c(1,1,1))
      D_xyz.ssv = matrix(rep(D_xyz.ssv,3),nRows,nCols-1)
      
      E <- acos(((xyz.coords.matrix/xyz.coords.ssv)*(D_xyz/D_xyz.ssv))%*%c(1,1,1))#kut izmedu normale tangencijalne ravnine i zrake prijemnik-satelit.
      ele <- apply(E,1,function.ele)
      
      Wii = 1/(sin(ele))^2 
      print("ele")
      print(ele)
      
    }else if(option==2){
      D_xyz = -P[1:nRows,1:(nCols-1)]/DD[1:nRows,1:(nCols-1)]
      D_xyz.ssv <- sqrt(D_xyz**2%*%c(1,1,1))
      D_xyz.ssv = matrix(rep(D_xyz.ssv,3),nRows,nCols-1)
      
      E <- acos(((xyz.coords.matrix/xyz.coords.ssv)*(D_xyz/D_xyz.ssv))%*%c(1,1,1))
      ele <- apply(E,1,function.ele)
      
      Wii = a_ele**2+b_ele**2/(sin(ele))^2
      print("ele")
      print(ele)
    }else if(option == 3){ 
      D_xyz = -P[1:nRows,1:(nCols-1)]/DD[1:nRows,1:(nCols-1)]
      D_xyz.ssv <- sqrt(D_xyz**2%*%c(1,1,1))
      D_xyz.ssv = matrix(rep(D_xyz.ssv,3),nRows,nCols-1)
      
      E <- acos(((xyz.coords.matrix/xyz.coords.ssv)*(D_xyz/D_xyz.ssv))%*%c(1,1,1))
      ele <- apply(E,1,function.ele)
      
      Wii = a**2/(sin(ele)+psi)^2
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
    } else {
      print("Nepoznata opcija rješavanja sustava jednadžbi.")
      exit(-1)
    }
    
    x_0 <- x_0 + dx*c(change,TRUE)
    delt <- dx*c(change,TRUE)
  }else{
    #delt <- -(1/2)*svd.inverse(PI)%*%PX
    delt <- qr.coef(qr(A_iter), b)
    x_0 = x_0 + delt #(x,y,z,c*dT)
  }
  
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

plot(iter, log10(abs(d.x)),type = 'l', col = 'red', main = c('LSA vrijeme izvršavanja u [s]=', round(timediff, digits = 2)), xlab = 'Broj iteracija', ylab = 'log10 delta_x po komponentama [m]')
lines(iter, log10(abs(d.y)), type = 'l', col = 'green')
lines(iter, log10(abs(d.z)), type = 'l', col = 'blue')
legend(
  "topright", legend=c("x","y", "z"), col=c("red", "green", "blue"),   pch=15
)

plot(iter, (abs(d.x)), type = 'l', col = 'red', main = c('LSA vrijeme izvršavanja u [s]=', round(timediff, digits = 2)), xlab = 'Broj iteracija', ylab = 'delta_x po komponentama [m]')
lines(iter, (abs(d.y)), type = 'l', col = 'green')
lines(iter, (abs(d.z)), type = 'l', col = 'blue')
legend(
  "topright", legend=c("x","y", "z"), col=c("red", "green", "blue"),   pch=15
)

iter <- err$V1
xx <- err$V2
yy <- err$V3
zz <- err$V4

plot(iter, log10(abs(xx)), 
     ylim = c(min(min(log10(abs(xx))),min(log10(abs(yy))),min(log10(abs(zz)))), max(max(log10(abs(xx))),max(log10(abs(yy))),max(log10(abs(zz))))),
     type = 'l', col = 'red', main = 'LSA odstupanje procjene položaja \n od stvarne vrijednosti', xlab = 'Broj iteracija', ylab = 'log10 pogrešaka [m]')
lines(iter, log10(abs(yy)), type = 'l', col = 'green')
lines(iter, log10(abs(zz)), type = 'l', col = 'blue')
legend(
"topright", legend=c("x","y", "z"), col=c("red", "green", "blue"),   pch=15
)

plot(iter, (abs(xx)), 
     ylim = c(min(min(abs(xx)),min(abs(yy)),min(abs(zz))), max(max(abs(xx)),max(abs(yy)),max(abs(zz)))),
     type = 'l', col = 'red', main = 'LSA odstupanje procjene položaja \n od stvarne vrijednosti', xlab = 'Broj iteracija', ylab = 'pogreška [m]')
lines(iter, (abs(yy)), type = 'l', col = 'green')
lines(iter, (abs(zz)), type = 'l', col = 'blue')
legend(
  "topright", legend=c("x","y", "z"), col=c("red", "green", "blue"),   pch=15
)

file.remove('razmakIteracija.txt','stvarnoOdstupanje.txt')