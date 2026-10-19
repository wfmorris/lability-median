# Make panels of Fig. 2 for manuscript entitled "Testing for demographic lability: 
# what should the baseline environment be?" 
# by William Morris and Daniel Doak
# "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

# Fig 2a plots the adult survival rate as power functions of driver z

# Fig 2b plots asymmetric Beta driver distributions with different amounts
# of right skew

# Fig 2c plots the combinations of curvature of the vital rate function and environmental
# driver variance (equivalently, skewness) at which the stochastic growth rate
# exceeds the growth rate that would result if the vital rate were fixed at its
# value at the mean or the median of the driver. Specifically, here we consider
# all vital rates to be fixed except for adult survival, 
# which follows a power function of driver z, with z following a right-skewed Beta
# distribution with mean 0.3 - thus greater Beta variance means greater right skew.
# Here we set the parameters of the adult survival function so that its value at
# mean(z) yields a deterministic population growth rate of 1 when the driver is 
# fixed at z=mean(z)

rm(list=ls(all=TRUE))
graphics.off()

# beta mean and variance of driver z
mu=0.3  # if mu<0.5, z is right-skewed
V.max = mu*(1-mu) # constraint on V for the beta

nvars=20  # number of values of the Beta variance to use
Vs=seq(0.01*V.max,0.2*V.max,0.19*V.max/(nvars-1)) # restrict variance to range where z distribution is unimodal

# use desired mean mu and variance V to get shape parameters for the beta 
as=mu*(mu*(1-mu)/Vs - 1)
bs=(1-mu)*(mu*(1-mu)/Vs - 1)

dz=0.002
z=seq(0,1,dz) # driver levels for plotting
nz=length(z) # number of driver levels for plotting

sj=0.1  # juvenile survival

# vital rate function for sa(z) is the power function sa(z)=z^B
Bs=seq(0.6,1.2,0.05)  # power function exponent; sa(z) concave if B<1, convex if B>1, linear if B=1
nBs=length(Bs)  # number of vital rate functions
sa.mu=mu^Bs     # adult survival at mean(z)
F=(1-sa.mu)/sj   # calculate fecundity values that yield a population growth rate 
                 # (i.e., dominant eigenvalue of the matrix [0 F;sj sa]) of 1 when 
                 # adult survival is at its value at z=mean(z), for each value of vital rate curvature

# verify that lambda=1 if z fixed at mean(z) (=mu)
lam.mean=array(0,nBs)
for (nn in 1:nBs) {
  lam.mean[nn]=max(eigen(matrix(c(0,sj,F[nn],sa.mu[nn]),2,2))$val)
}

# Fig. 2A

# plot adult survival functions
sa.z=matrix(0,nz,nBs)
for(col in 1:nBs) sa.z[,col]=z^Bs[col]
windows()
matplot(z,sa.z,type='l',
        xlab='Driver, z',ylab='Adult survival',lwd=3,cex.axis=2.5,lty='solid',xaxt='n')
lines(c(.3,.3),c(0,1),lwd=4,lty='dashed')

# Fig. 2B

# plot the Beta distributions of z
beta.pdf=matrix(0,nz,nvars)
for(col in 1:nvars) beta.pdf[,col]=dbeta(z,as[col],bs[col])
windows()
matplot(z,beta.pdf,type='l',
        xlab='Driver, z',ylab='Probability density',lwd=3,cex.axis=2.5,lty='solid')
lines(c(.3,.3),c(0,9),lwd=4,lty='dashed')

Ms=qbeta(0.5,as,bs)  # compute median z for all driver distributions

# matrices to store results; rows = sa functions, columns = z distributions
sa.M=lambda.M=lambda.s=matrix(0,nBs,nvars)

# compute adult survival at each MEDIAN z
for(row in 1:nBs) sa.M[row,]=Ms^Bs[row]

# compute deterministic lambda when z is fixed at its MEDIAN, for all sa(z) functions
for(row in 1:nBs){
  for(col in 1:nvars){
    P=matrix(c(0,sj,F[row],sa.M[row,col]),2,2)  # projection matrix at MEDIAN z
    lambda.M[row,col]=eigen(P)$val[1]
  }
}

# compute stochastic lambda when z varies, for all z distributions and sa(z) functions
tmax=5E4 # years to compute stochastic lambda
buffer=1E3  # number of initial years to discard before convergence to stationary distribution is reached

nreps=100 # number of times to compute stochastic lambda for each distribution of z
Bcrit=matrix(0,nreps,nvars) # array to store value of power function exponent where lambda.s first exceeds lambda.M
flag=array(0,nvars) # array to store whether lambda.s exceeds lambda.M even at lowest B explored

for(rep in 1:nreps){
  print(rep)
  for(col in 1:nvars){ # loop over distributions of driver z
    zs=rbeta(tmax,as[col],bs[col]) # generate tmax random values of z; use with each sa function
    for(row in 1:nBs){  # loop over curvature values of the adult survival function
      lams=array(0,tmax) # array to store annual population growth rates
      n0=matrix(c(0.5,0.5),2,1)  # starting population vector - uniform stage structure
      sa.t=zs^Bs[row] # adult survival probabilities for each year
      for(t in 1:tmax){
        P=matrix(c(0,sj,F[row],sa.t[t]),2,2)  # annual projection matrix
        n1=P %*% n0
        lams[t]=sum(n1)  
        n0=n1/sum(n1) # renormalize population vector each year
      }
      lambda.s[row,col]=exp(mean(log(lams[buffer:tmax]))) # compute geometric mean growth rate
    }
  }
  delta.lambda = lambda.s-lambda.M # difference between stochastic lambda and lambda at median z
  for(col in 1:nvars){
    index=max(which(delta.lambda[,col]<0))  
    if(index == - Inf) {
      flag[col]=1
      index = 1  # if delta.lambda doesn't change sign, save the lowest B used
    }
    Bcrit[rep,col]=Bs[index] # find lowest B where stochastic lambda exceeds lambda.M
  }
}
meanBcrit=colMeans(Bcrit)

# Fig. 2C

windows()
plot(Vs,meanBcrit,type='b',pch=NA,
     xlim=c(Vs[1],Vs[nvars]),ylim=c(0.4,1.2),xlab="Driver variance",ylab="Annual survival exponent b",
     lwd=3,cex.axis=2.5,lty='solid')
for(s in 2:length(meanBcrit)) lines(c(Vs[s-1],Vs[s]),c(meanBcrit[s-1],meanBcrit[s]),lty='solid',lwd=3)
lines(c(Vs[1],Vs[nvars]),c(1,1),lty='dashed',lwd=3)