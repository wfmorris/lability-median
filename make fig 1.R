# Make panels of Fig 1 for manuscript entitled "Testing for demographic lability: 
# what should the baseline environment be?" by William Morris and Daniel Doak
# "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

# Fig 1 plots nine combinations of curvature of a vital rate function (convex, 
# linear, and concave) and skewness of the environmental driver distribution 
# (left-skewed, symmetric, and right-skewed). Also plotted are the driver median 
# (as a dark blue vertical dashed line) and mean (as a cyan vertical dashed line). 
# The script also computes the proportional change in the arithmetic mean of the 
# vital rate when the driver varies according to the driver distribution shown
# relative to the fixed vital rate value at the driver median (L.med)
# and at the driver mean (L.mean). L is referred to as the "lability index" in 
# the main text


graphics.off()

Mu.mean=Mu.med=Mu.v=array(0,9)  # Mu.mean stores adult survival at mean z
                                # Mu.med stores adult survival at median z
                                # Mu.v stores arithmetic mean adult survival when z varies

# ***** Panel a *****

# mean and variance of driver z
M=0.7  # left-skewed
V.max = M*(1-M) # constraint on V for the Beta
V=0.03

# use desired mean M and variance V to get shape parameters for the Beta 
a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)
z.med=qbeta(0.5,a,b)  # get median z

z=seq(0,1,.002) # driver levels

B=2  # convex vital rate/driver relationship (VRDR)

sa.z=z^B  # vital rate function

drivers=rbeta(1E8,a,b) # sample from driver distribution
Mu.mean[1] = M^B       # vital rate at mean of driver
Mu.med[1] = z.med^B       # vital rate at median of driver
Mu.v[1] = mean(drivers^B ) # mean vital rate when driver varies

windows()
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(z, dbeta(z,a,b),type='l',xlab='',ylab='',lwd=3,lty='solid',las=1,cex.axis=2,xaxt='n') # first plot
axis(1,xlim=c(0,1),labels=FALSE,cex.axis=2,las=1)
axis(2,ylim=c(0,2.5),labels=TRUE,cex.axis=2,las=1)
par(new = TRUE)
plot(z, sa.z, type = "l", bty = "n", xlab = "", ylab = "",lwd=3,cex.axis=2,lty='solid',col='red',axes=FALSE)
axis(side=4, at = pretty(range(z)),col.axis='red',labels=FALSE,cex.axis=2,las=1)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3,col='blue')
lines(c(M,M),c(0,2.5),lty='dashed',lwd=3,col='cyan')

# ***** Panel b *****

# mean and variance of driver z
M=0.5  # symmetric
V.max = M*(1-M) # constraint on V for the Beta
V=0.03

# use desired mean M and variance V to get shape parameters for the Beta 
a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)
z.med=qbeta(0.5,a,b)  # get median z

z=seq(0,1,.002) # driver levels

B=2  # convex VRDR

sa.z=z^B 

drivers=rbeta(1E8,a,b) # sample from driver distribution
Mu.mean[2] = M^B       # vital rate at mean of driver
Mu.med[2] = z.med^B          # vital rate at median of driver
Mu.v[2] = mean(drivers^B ) # mean vital rate when driver varies

windows()
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(z, dbeta(z,a,b),type='l',xlab='',ylab='',lwd=3,lty='solid',las=1,cex.axis=2,xaxt='n',yaxt='n') # first plot
axis(1,xlim=c(0,1),labels=FALSE,cex.axis=2,las=1)
axis(2,ylim=c(0,2.5),labels=FALSE,cex.axis=2,las=1)
par(new = TRUE)
plot(z, sa.z, type = "l", bty = "n", xlab = "", ylab = "",lwd=3,cex.axis=2,lty='solid',col='red',axes=FALSE)
axis(side=4, at = pretty(range(z)),col.axis='red',labels=FALSE,cex.axis=2,las=1)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3,col='blue')
lines(c(M,M),c(0,2.5),lty='dashed',lwd=3,col='cyan')

# ***** Panel c *****

# mean and variance of driver z
M=0.3  # right-skewed
V.max = M*(1-M) # constraint on V for the beta
V=0.03

# use desired mean M and variance V to get shape parameters for the Beta 
a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)
z.med=qbeta(0.5,a,b)  # get median z

z=seq(0,1,.002) # driver levels

B=2  # convex VRDR

sa.z=z^B 

drivers=rbeta(1E8,a,b) # sample from driver distribution
Mu.mean[3]= M^B       # vital rate at mean of driver
Mu.med[3] = z.med^B       # vital rate at median of driver
Mu.v[3] = mean(drivers^B ) # mean vital rate when driver varies

windows()
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(z, dbeta(z,a,b),type='l',xlab='',ylab='',lwd=3,lty='solid',las=1,cex.axis=2,xaxt='n',yaxt='n') # first plot
axis(1,xlim=c(0,1),labels=FALSE,cex.axis=2,las=1)
axis(2,ylim=c(0,2.5),labels=FALSE,cex.axis=2,las=1)
par(new = TRUE)
plot(z, sa.z, type = "l", bty = "n", xlab = "", ylab = "",lwd=3,cex.axis=2,lty='solid',col='red',axes=FALSE)
axis(side=4, at = pretty(range(z)),col.axis='red',labels=TRUE,cex.axis=2,las=1)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3,col='blue')
lines(c(M,M),c(0,2.5),lty='dashed',lwd=3,col='cyan')

# ***** Panel d *****

# mean and variance of driver z
M=0.7  # left-skewed
V.max = M*(1-M) # constraint on V for the beta
V=0.03

# use desired mean M and variance V to get shape parameters for the Beta 
a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)
z.med=qbeta(0.5,a,b)   # get median z

z=seq(0,1,.002) # driver levels

B=1  # linear VRDR

sa.z=z^B 

drivers=rbeta(1E8,a,b) # sample from driver distribution
Mu.mean[4] = M^B       # vital rate at mean of driver
Mu.med[4] = z.med^B       # vital rate at median of driver
Mu.v[4] = mean(drivers^B ) # mean vital rate when driver varies

windows()
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(z, dbeta(z,a,b),type='l',xlab='',ylab='',lwd=3,lty='solid',las=1,cex.axis=2,xaxt='n',yaxt='n') # first plot
axis(1,xlim=c(0,1),labels=FALSE,cex.axis=2,las=1)
axis(2,ylim=c(0,2.5),labels=TRUE,cex.axis=2,las=1)
par(new = TRUE)
plot(z, sa.z, type = "l", bty = "n", xlab = "", ylab = "",lwd=3,cex.axis=2,lty='solid',col='red',axes=FALSE)
axis(side=4, at = pretty(range(z)),col.axis='red',labels=FALSE,cex.axis=2,las=1)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3,col='blue')
lines(c(M,M),c(0,2.5),lty='dashed',lwd=3,col='cyan')

# ***** Panel e *****

# beta mean and variance of driver z
M=0.5  # symmetric
V.max = M*(1-M) # constraint on V for the beta
V=0.03

# use desired mean M and variance V to get shape parameters for the Beta 
a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)
z.med=qbeta(0.5,a,b)   # get median z

z=seq(0,1,.002) # driver levels

B=1  # linear VRDR

sa.z=z^B 

drivers=rbeta(1E8,a,b) # sample from driver distribution
Mu.mean[5]  =M^B       # vital rate at mean of driver
Mu.med[5] = z.med^B          # vital rate at median of driver
Mu.v[5] = mean(drivers^B ) # mean vital rate when driver varies

windows()
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(z, dbeta(z,a,b),type='l',xlab='',ylab='',lwd=3,lty='solid',las=1,cex.axis=2,xaxt='n',yaxt='n') # first plot
axis(1,xlim=c(0,1),labels=FALSE,cex.axis=2,las=1)
axis(2,ylim=c(0,2.5),labels=FALSE,cex.axis=2,las=1)
par(new = TRUE)
plot(z, sa.z, type = "l", bty = "n", xlab = "", ylab = "",lwd=3,cex.axis=2,lty='solid',col='red',axes=FALSE)
axis(side=4, at = pretty(range(z)),col.axis='red',labels=FALSE,cex.axis=2,las=1)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3,col='blue')
lines(c(M,M),c(0,2.5),lty='dashed',lwd=3,col='cyan')

# ***** Panel f *****

# mean and variance of driver z
M=0.3  # right-skewed
V.max = M*(1-M) # constraint on V for the beta
V=0.03

# use desired mean M and variance V to get shape parameters for the Beta 
a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)
z.med=qbeta(0.5,a,b)  # get median z

z=seq(0,1,.002) # driver levels

B=1  # linear VRDR

sa.z=z^B 

drivers=rbeta(1E8,a,b) # sample from driver distribution
Mu.mean[6]=M^B       # vital rate at mean of driver
Mu.med[6] = z.med^B          # vital rate at median of driver
Mu.v[6] = mean(drivers^B ) # mean vital rate when driver varies

windows()
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(z, dbeta(z,a,b),type='l',xlab='',ylab='',lwd=3,lty='solid',las=1,cex.axis=2,xaxt='n',yaxt='n') # first plot
axis(1,xlim=c(0,1),labels=FALSE,cex.axis=2,las=1)
axis(2,ylim=c(0,2.5),labels=FALSE,cex.axis=2,las=1)
par(new = TRUE)
plot(z, sa.z, type = "l", bty = "n", xlab = "", ylab = "",lwd=3,cex.axis=2,lty='solid',col='red',axes=FALSE)
axis(side=4, at = pretty(range(z)),col.axis='red',labels=TRUE,cex.axis=2,las=1)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3,col='blue')
lines(c(M,M),c(0,2.5),lty='dashed',lwd=3,col='cyan')

# ***** Panel g *****

# mean and variance of driver z
M=0.7  # left-skewed
V.max = M*(1-M) # constraint on V for the beta
V=0.03

# use desired mean M and variance V to get shape parameters for the Beta 
a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)
z.med=qbeta(0.5,a,b)  # get median z

z=seq(0,1,.002) # driver levels

B=0.7  # concave VRDR

sa.z=z^B 

drivers=rbeta(1E8,a,b) # sample from driver distribution
Mu.mean[7]=M^B       # vital rate at mean of driver
Mu.med[7] = z.med^B          # vital rate at median of driver
Mu.v[7] = mean(drivers^B ) # mean vital rate when driver varies

windows()
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(z, dbeta(z,a,b),type='l',xlab='',ylab='',lwd=3,lty='solid',las=1,cex.axis=2,xaxt='n',yaxt='n') # first plot
axis(1,xlim=c(0,1),labels=TRUE,cex.axis=2,las=1)
axis(2,ylim=c(0,2.5),labels=TRUE,cex.axis=2,las=1)
par(new = TRUE)
plot(z, sa.z, type = "l", bty = "n", xlab = "", ylab = "",lwd=3,cex.axis=2,lty='solid',col='red',axes=FALSE)
axis(side=4, at = pretty(range(z)),col.axis='red',labels=FALSE,cex.axis=2,las=1)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3,col='blue')
lines(c(M,M),c(0,2.5),lty='dashed',lwd=3,col='cyan')

# ***** Panel h *****

# mean and variance of driver z
M=0.5  # symmetric
V.max = M*(1-M) # constraint on V for the beta
V=0.03

# use desired mean M and variance V to get shape parameters for the Beta 
a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)
z.med=qbeta(0.5,a,b)  # get median z

z=seq(0,1,.002) # driver levels

B=0.7  # concave VRDR

sa.z=z^B 

drivers=rbeta(1E8,a,b) # sample from driver distribution
Mu.mean[8] = M^B       # vital rate at mean of driver
Mu.med[8] = z.med^B       # vital rate at median of driver
Mu.v[8] = mean(drivers^B ) # mean vital rate when driver varies

windows()
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(z, dbeta(z,a,b),type='l',xlab='',ylab='',lwd=3,lty='solid',las=1,cex.axis=2,xaxt='n',yaxt='n') # first plot
axis(1,xlim=c(0,1),labels=TRUE,cex.axis=2,las=1)
axis(2,ylim=c(0,2.5),labels=FALSE,cex.axis=2,las=1)
par(new = TRUE)
plot(z, sa.z, type = "l", bty = "n", xlab = "", ylab = "",lwd=3,cex.axis=2,lty='solid',col='red',axes=FALSE)
axis(side=4, at = pretty(range(z)),col.axis='red',labels=FALSE,cex.axis=2,las=1)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3,col='blue')
lines(c(M,M),c(0,2.5),lty='dashed',lwd=3,col='cyan')

# ***** Panel i *****

# mean and variance of driver z
M=0.3  # right-skewed
V.max = M*(1-M) # constraint on V for the beta
V=0.03

# use desired mean M and variance V to get shape parameters for the Beta 
a=M*(M*(1-M)/V - 1)
b=(1-M)*(M*(1-M)/V - 1)
z.med=qbeta(0.5,a,b)  # get median z

z=seq(0,1,.002) # driver levels

B=0.7  # concave VRDR

sa.z=z^B 

drivers=rbeta(1E8,a,b) # sample from driver distribution
Mu.mean[9]=M^B       # vital rate at mean of driver
Mu.med[9] = z.med^B          # vital rate at median of driver
Mu.v[9] = mean(drivers^B ) # mean vital rate when driver varies

windows()
par(mar = c(5, 4, 4, 4) + 0.5)  # Leave space for z axis
plot(z, dbeta(z,a,b),type='l',xlab='',ylab='',lwd=3,lty='solid',las=1,cex.axis=2,xaxt='n',yaxt='n') # first plot
axis(1,xlim=c(0,1),labels=TRUE,cex.axis=2,las=1)
axis(2,ylim=c(0,2.5),labels=FALSE,cex.axis=2,las=1)
par(new = TRUE)
plot(z, sa.z, type = "l", bty = "n", xlab = "", ylab = "",lwd=3,cex.axis=2,lty='solid',col='red',axes=FALSE)
axis(side=4, at = pretty(range(z)),col.axis='red',labels=TRUE,cex.axis=2,las=1)
lines(c(z.med,z.med),c(0,2.5),lty='dashed',lwd=3,col='blue')
lines(c(M,M),c(0,2.5),lty='dashed',lwd=3,col='cyan')

# compute lability indices for all panels, relative to the driver mean and median
L.mean=(Mu.v-Mu.mean)/Mu.mean # lability index using the mean as baseline
L.med=(Mu.v-Mu.med)/Mu.med # lability index using the median as baseline
