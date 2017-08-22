x = which(min(abs(xx)) == abs(xx))
y = which(min(abs(yy)) == abs(yy))
z = which(min(abs(zz)) == abs(zz))

xval = abs(xx)[which(min(abs(xx)) == abs(xx))]
yval = abs(yy)[which(min(abs(yy)) == abs(yy))]
zval = abs(zz)[which(min(abs(zz)) == abs(zz))]
cat(c(x,y,z,xval,yval,zval),' \r',file="x_y_z_xval_yval_zval.txt", append=TRUE)