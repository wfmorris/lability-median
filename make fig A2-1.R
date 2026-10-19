# Make Fig. A2.1 for manuscript entitled "Testing for demographic lability: 
# what should the baseline environment be?" 
# by William Morris and Daniel Doak
# "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

# The figure plot adult survival vs driver z using the modified Beta probability 
# density function given my eq. (A2.1) in appendix 2 of the manuscript

graphics.off()

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "red", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
alphas=c(2,4,8,2,4,8,4,8,16)
betas=c(2,4,8,1,1,1,2,4,8)
# to plot just one function: 
# alphas=c(8); betas=c(1)
na=length(alphas)
modes=(alphas-1)/(alphas+betas-2)
z=seq(0,1,.001)
nz=length(z)
y=matrix(0,nz,na)
for(c in 1:na){
  y[,c]=0.99*dbeta(z,alphas[c],betas[c])/dbeta(modes[c],alphas[c],betas[c])
}
windows()
par(mar = c(5, 4, 4, 4) + 1)  # Leave space for z axis
matplot(z,y,type='l',lty='solid',lwd=3,cex.axis=2,xlab='Environment, z',ylab='Adult survival',
        cex.lab=2,col=safe_colorblind_palette)
legend("topleft",c('2,2','4,4','8,8','2,1','4,1','8,1','4,2','8,4','16,8'),lty='solid',lwd=3,
       col=safe_colorblind_palette)