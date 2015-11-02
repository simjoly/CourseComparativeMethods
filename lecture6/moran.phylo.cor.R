#-------------------------------------------------
#
# Moran's I phylogenetic autocorrelograms
#
# Simon Joly, Montreal Botanical Garden, 2008-2012
#
#-------------------------------------------------

moran.phylo.cor <- function(XX.vec,YY.mat,breaks=4,p.value=0.05,plot=c("FALSE","TRUE"))

# The function moran.phylo.cor performs and plot a Moran's I phylogenetic autocorrelogram. The method used for the autocorrelogram follows that described in chapter 13 of Legendre and Legendre (1998) for Mantel spatial correlogram.

# The function requires the package "ape" to be installed for the calculation of the Moran I statistic
#
# INPUT Parameters:
#
# XX.vec:        vector of traits
# YY.mat:        Phylogenetic distance matrix
# breaks:        number of breaks in the data
# p.value:       The cutoff p-value used for statistical significance 
# plot:          A string indicating if the Moran's I correlogram should be
#                plotted (default) or not.
#
#
# The function returns an OUTPUT list containing the following ELEMENTS:
#
# values:         The mantel correlation values for each distance class
# significance:   The P values for mantel tests in each class
# sd:        	    Standard deviation
# breaks:		      the size of the intervals
# mids:			      the median values of each interval

{
require(ape)

breaks = breaks[1]
intervals = (max(YY.mat) - min(YY.mat)) / breaks
the.breaks = min(YY.mat)
for (i in 1:breaks) {
	the.breaks <- cbind(the.breaks,min(YY.mat)+(intervals*i))
}

# Create vectors to store the mantel values and their significance
moran.values <- rep(NA, (length(the.breaks)-1))
moran.significance <- rep(1, (length(the.breaks)-1))
moran.sd <- rep(1, (length(the.breaks)-1))

a=length(the.breaks)-1
i=1

while (a){
	YY.mat1 <- YY.mat
    z=nrow(YY.mat)
    x=1
    while (z){
        zz=nrow(YY.mat)
        y=1
        while (zz){
            if ( (YY.mat[x,y] <= the.breaks[i]) | (YY.mat[x,y] > the.breaks[(i+1)])) YY.mat1[x,y]=0
            else YY.mat1[x,y]=1
            y=y+1
            zz=zz-1
        }
        x=x+1
        z=z-1
    }
	# skip if no distances in the interval
   	if (max(YY.mat1)==0) {
   		i=i+1
   		a=a-1
   		next
   	}
    moran.result <- Moran.I(XX.vec,YY.mat1)
  moran.values[i] <- moran.result$observed
	moran.significance[i] <- moran.result$p.value
    moran.sd[i] <- moran.result$sd
    i=i+1
    a=a-1
}
mids <- numeric(breaks)
for (i in 1:breaks) {
	mids[i]=the.breaks[i]+intervals/2
}
if(plot!=FALSE){
	value.max <- max(moran.values[!is.na(moran.values)])
	value.min <- min(moran.values[!is.na(moran.values)])
	xaxis.min = min(the.breaks)
	xaxis.max = max(the.breaks)
	#i=1
	op <- par(mar=c(2, 4, 4, 2),mgp=c(2.5,1,0))
	plot(moran.values,type="n",bty="n",xlim=c(xaxis.min,xaxis.max),ylim=c(value.min,value.max),
		cex.axis=0.8,xaxt="n",font.lab=2,main="Moran's I Phylogenetic correlogram",
		xlab="",ylab="Moran's I")
	for (i in 1:length(moran.values)){
	  if (is.na(moran.values[i])) next
	  if (moran.significance[i] <= p.value) points(mids[i],moran.values[i],pch=22,cex=2,bg="black")
	  else points(mids[i],moran.values[i],pch=22,cex=2,bg="white")
	  #i=i+1
	  }
	abline(h=0,lty=2)
	par(op)
}

# Return the results
return(list(values=moran.values, significance=moran.significance, sd=moran.sd,breaks=as.vector(the.breaks),mids=mids))
}
