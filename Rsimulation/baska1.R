# S1
x1 = 7766188.44
y1 = -21960535.34
z1 = 12522838.56
# S2
x2 = -25922679.66
y2 = -6629461.28
z2 = 31864.37
# S3
x3 = -5743774.02
y3 = -25828319.92
z3 = 1692757.72
# S4
x4 = -2786005.69
y4 = -15900725.80
z4 = 21302003.49
# Pseudoranges
rho1 = 22228206.42
rho2 = 24096139.11
rho3 = 21729070.63
rho4 = 21259581.09

########################
x1 = 2088202.299
y1 = -11757191.370
z1 = 25391471.881
# S2
x2 = 11092568.240
y2 = -14198201.090
z2 = 21471165.950
# S3
x3 = 35606984.591
y3 = -4447027.239
z3 = 9101378.572
# S4
x4 = 3966929.048
y4 = 7362851.831
z4 = 26388447.172
# Pseudoranges
rho1 = 23204698.51 #calculated from PRN codes
rho2 = 21585835.37
rho3 = 31364260.01
rho4 = 24966798.73


# TO DO: Reading sat co-ordinates and pseudoranges from files using:
# data <- read.csv("G:\\satellites.dat", header = FALSE, sep = ' ')

# Physical constants
c = 299792458 # 300000000 #na prvom linku je 1
R = 6378137

# Initial values of errors
delta_x = 11
delta_y = 11
delta_z = 11
delta_t = 11

# Initial values for iteration
xk = 0
yk = 0
zk = 0
tk = 0

#
C = diag(4)
delta = 3;

# Loop until error is less than or equal to 1 m - provisional
while (abs(delta) > 1) { # TO DO: Assure that all the errors are below 1 m
  #simple iteration as in https://www.u-blox.com/sites/default/files/products/documents/GPS-Compendium_Book_%28GPS-X-02007%29.pdf
  #expected error is 0.
  R1 = sqrt((xk - x1) ** 2 + (yk - y1) ** 2 + (zk - z1) ** 2)
  R2 = sqrt((xk - x2) ** 2 + (yk - y2) ** 2 + (zk - z2) ** 2)
  R3 = sqrt((xk - x3) ** 2 + (yk - y3) ** 2 + (zk - z3) ** 2)
  R4 = sqrt((xk - x4) ** 2 + (yk - y4) ** 2 + (zk - z4) ** 2)
  
  ax1 = (xk - x1) / R1 #U prvom link rho_j i R_j su zamjenjeni
  ay1 = (yk - y1) / R1
  az1 = (zk - z1) / R1
  
  ax2 = (xk - x2) / R2
  ay2 = (yk - y2) / R2
  az2 = (zk - z2) / R2
  
  ax3 = (xk - x3) / R3
  ay3 = (yk - y3) / R3
  az3 = (zk - z3) / R3
  
  ax4 = (xk - x4) / R4
  ay4 = (yk - y4) / R4
  az4 = (zk - z4) / R4
  
  delta_rho = matrix(c(rho1 - R1, rho2 - R2,
                       rho3 - R3, rho4 - R4),
                     nrow = 4, ncol = 1, byrow = TRUE)
  # Matrix H (or G) 
  #TO DO: Least-Sqaures Adjustment for n (number of satellites) > 4 in a form:
  #x = (GT G)^(-1) GT y
  H = matrix(c(ax1, ay1, az1, c,
               ax2, ay2, az2, c,
               ax3, ay3, az3, c,
               ax4, ay4, az4, c),
              nrow = 4, ncol = 4, byrow = TRUE)
  
  # Inverse of matrix H
  Hinv = solve(H)
  
  # Dot product of inversed matrix H and calculated pseudoranges
  dx = Hinv %*% C %*% delta_rho
  
  # Errors
  delta_x = dx[1]
  delta_y = dx[2]
  delta_z = dx[3]
  delta_t = dx[4]
  
  # Pseudorange calculation
  rho1 = R1 + ax1 * delta_x + ay1 * delta_y + az1 * delta_z
  rho2 = R2 + ax2 * delta_x + ay2 * delta_y + az2 * delta_z
  rho3 = R3 + ax3 * delta_x + ay3 * delta_y + az3 * delta_z
  rho4 = R4 + ax4 * delta_x + ay4 * delta_y + az4 * delta_z
  
  # Calculating coordinates for next iteration step
  xk = xk + delta_x #formula 6.7 iz prvog linka
  yk = yk + delta_y
  zk = zk + delta_z
  tk = tk + delta_t
  delta = max(c(delta_x,delta_y,delta_z,delta_t));  

  print ("This iteration position estimation: ")
  print (c(xk, yk, zk, tk))
  print("Each coordinate error: ");
  print (c(delta_x, delta_y, delta_z))
}

# Calculating latitude and longitude and converting from radians to degrees
lat = asin(zk / R) * (180 / pi)
lon = atan2(yk, xk) * (180 / pi)

print (c(lat, lon))

#TO DO: data collection for graphical presentation of:
# - efficiency of the iteration algorithm
# - position estimation accuracy

