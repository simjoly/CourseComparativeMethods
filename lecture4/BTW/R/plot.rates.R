plot.rates = function (model, estimates=T, type="n", xaxt="n", yaxt="n", xlab="", ylab="", ...) {
	
########################
### CHECK MODEL TYPE ###
########################

	if (class(model) %in% c("DiscreteML", "DiscreteMCMC") == F) {
		stop("'plot.rates' requires a Discrete BayesTraits object")
	}
	
############
### PLOT ###
############

mat = model$RateMatrix[[1]]
rates <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),2) -> labs
if (estimates == F) {labs=c("q12", "q21", "q24", "q42", "q43", "q34", "q31", "q13")}
plot(c(0,100), c(0,100), type=type, xaxt=xaxt, yaxt=yaxt, xlab=xlab, ylab=ylab, ...)
text(x=c(20, 80, 20, 80), y=c(80, 80, 20, 20), labels=c("00", "01", "10", "11"), cex=4)
text(x=c(50, 50, 93, 67, 50, 50, 7, 33), y=c(93, 67, 50, 50, 7,33, 50, 50), labels=labs, cex=0.75)
arrows(x0=c(35, 65, 85, 75, 65, 35, 15, 25), y0=c(85, 75, 65, 35, 15, 25, 35, 65), x1=c(65, 35, 85, 75, 35, 65, 15, 25), y1=c(85, 75, 35, 65, 15, 25, 65, 35), lwd=rates/max(rates, na.rm=T)*15)

}



