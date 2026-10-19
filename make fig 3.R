# Make panels of Fig. 3 for manuscript entitled "Testing for demographic lability: 
# what should the baseline environment be?" 
# by William Morris and Daniel Doak
# "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

# Fig 3a plots the stochastic population growth rate lambda.s - the measure of
# fitness - against the centering parameter of the adult survival function given 
# by eq. (A2.1) in appendix 2 of the manuscript.

# Fig. 3b plots the initial adult survival function, the adult survival functions
# centered on the driver mean and median, and the optimum survival rate 
# function that maximizes lambda.s

rm(list=ls(all=TRUE))
graphics.off()

tmax=5e5 # "years" for computing stochastic lambda

sj=0.1  # fixed juvenile survival for the projecton matrix in eq. (2) in the main text

# ******************************************************
# Here we consider a negatively skewed driver distribution, with mean mu.e > 0.5

mu.e = 0.7  # mean of the Beta environmental distribution
var.e = 0.02 # variance of the Beta environmental distribution
ae=mu.e*(mu.e*(1-mu.e)/var.e - 1) # alpha for the distribution of z
be=(1-mu.e)*ae/mu.e               # beta for the distribution of z
zs=rbeta(tmax,ae,be)   # random sample of environmental values z to use for all adult survival functions
med.e = qbeta(0.5,ae,be) # median of the Beta environmental distribution
z=seq(0,1,.001)   # range of driver z
p.e=dbeta(z,ae,be)
p.e=p.e/max(p.e)  #rescale distn. of z to its max, to compare to survival functions

# Initial adult survival function is convex with 'mean' a.0+b.0=0.8
a.0=4 
b.0=1
mode=(a.0-1)/(a.0+b.0-2)
scale.0=0.99/dbeta(mode,a.0,b.0)  # scaling parameter for the initial adult survival function
s.0=scale.0*dbeta(z,a.0,b.0)  # initial adult survival function 
Vs=a.0*b.0/((a.0+b.0)^2 * (a.0+b.0+1)) # "width" of the initial adult survival 
                                      # function = environmental tolerance, which
                                      # is assumed to be fixed

# Survival function centered on the mean z (mu.e) with the same "width" (var) as the initial function 
M=mu.e
a=M*(M*(1-M)/Vs - 1)
b=(1-M)*a/M
mode=(a-1)/(a+b-2)
scale=0.99/dbeta(mode,a,b)
s.mu=scale*dbeta(z,a,b)  # survival function centered on environmental mean

# Survival function centered on the median z with the same "width" (var) as the initial function 
M=med.e
a=M*(M*(1-M)/Vs - 1)
b=(1-M)*a/M
mode=(a-1)/(a+b-2)
scale=0.99/dbeta(mode,a,b)
s.med=scale*dbeta(z,a,b)  # survival function centered on environmental median

# First find value of fecundity F that yields lam.s = 1 for the initial convex 
# adult survival function with alpha = 4 and beta = 1

lam=array(0,tmax)

F=6.320 # desired value of F
A=matrix(c(0,sj,F,0),2,2) # constant parts of projection matrix
n1=array(0.5,2)  # starting population vector
for(t in 1:tmax){
  A[2,2]=scale.0*dbeta(zs[t],a.0,b.0) # get adult survival for the current z
  n2=A %*% n1
  lam[t]=sum(n2)
  n1=n2/sum(n2) # renormalize population vector
}
lam.s.0=exp(mean(log(lam))) # stochastic lambda is geometric mean of annual lambdas
print(lam.s.0) # verify that F yields stochastic lambda = 1

# Using F and sj, compute for the initial adult survival function the 
# deterministic lambdas at mean(z)=mu.e and median(z))=med.e
A[2,2]=scale.0*dbeta(mu.e,a.0,b.0) # get adult survival at mean(z)
lam.mu.e=eigen(A)$val[1]
A[2,2]=scale.0*dbeta(med.e,a.0,b.0) # get adult survival at median(z)
lam.med.e=eigen(A)$val[1]

# Fix "width" (variance) and only allow "center" (mean) of survival function to evolve
Ms=seq(.05,.95,.01)
as=Ms*(Ms*(1-Ms)/Vs - 1)
bs=(1-Ms)*as/Ms
I=which(round(as,digits=5)>1 & round(bs,digits=5)>=1) # find survival functions with alpha>1 and beta>=1, to avoid Inf survival
Ms=Ms[I]
nms=length(Ms)  # number of values of the centering parameter to consider
lam.s=array(0,nms)

for(m in 1:nms){
  print(m)
  Mm = Ms[m]
  a=round(Mm*(Mm*(1-Mm)/Vs - 1),digits=10)
  b=round((1-Mm)*a/Mm,digits=10)
  mode=(a-1)/(a+b-2)
  scale=0.99/dbeta(mode,a,b)
  A=matrix(c(0,sj,F,0),2,2) # constant parts of projection matrix
  n1=array(0.5,2)  # starting population vector
  for(t in 1:tmax){
    A[2,2]=scale*dbeta(zs[t],a,b) # get adult survival for the current z
    n2=A %*% n1
    lam[t]=sum(n2)
    n1=n2/sum(n2) # renormalize population vector
  }
  lam.s[m]=exp(mean(log(lam))) # stochastic lambda is geometric mean of annual lambdas
}
I=which(lam.s==max(lam.s))
lam.s.best=lam.s[I]
M.best=Ms[I]
ab=M.best*(M.best*(1-M.best)/Vs - 1)
bb=(1-M.best)*ab/M.best
mode=(ab-1)/(ab+bb-2)
scale=0.99/dbeta(mode,ab,bb)
s.best=scale*dbeta(z,ab,bb)

A[2,2]=scale*dbeta(mu.e,ab,bb) # get adult survival at mean(z)
lam.best.mu.e=eigen(A)$val[1] # deterministic lambda at mean(z) for opt surv curve
A[2,2]=scale*dbeta(med.e,ab,bb) # get adult survival at median(z)
lam.best.med.e=eigen(A)$val[1] # deterministic lambda at median(z) for opt surv curve

# plot the fitness function for the centering parameter z'
windows()
plot(Ms,lam.s,type='l',lwd=3,ylim=c(0,1.3),xlim=c(min(Ms),max(Ms)),
     xlab='Mean parameter of adult survival rate function',
     ylab='Stochastic lambda',cex.axis=2)
lines(mu.e*c(1,1),c(0,1.3),lty='dotted',lwd=3,col='blue')
lines(med.e*c(1,1),c(0,1.3),lty='dashed',lwd=3,col='blue')
lines(M.best*c(1,1),c(0,1.3),lty='solid',lwd=3,col='blue')
lines(c(min(Ms),max(Ms)),c(1,1),lty='dotted',lwd=2)

# compare 4 adult survival functions (initial and centered on mean(z), median(z), and optimal z) to the scaled environmental distribution
windows()
matplot(z,cbind(p.e,s.0,s.best,s.mu,s.med),type='l',lwd=3,col=c('black','red','blue','blue','blue'),
        lty=c('solid','solid','solid','dotted','dashed'),
        xlab='Environment, z',ylab='Survival rate or scaled density of z',
        cex.axis=2)