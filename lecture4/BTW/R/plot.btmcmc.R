# arg may be a vector with MCMC output or a BayesTraits MCMC object
# plot.singles only applies to vectors
# if arg is a vector, plots won't include the name of the parameter

plot.btmcmc = function (arg, plot.singles = FALSE) {

###################################
### CHECK FOR REQUIRED PACKAGES ###
###################################

	if (("coda" %in% rownames(installed.packages())) == F) {
		stop("You must install the package 'coda'.")
	}
	library(coda)	

#########################################
### if arg is a vector of MCMC output ###
#########################################
	
if (is.vector(arg) == TRUE && class(arg) %in% c("integer", "numeric")) {
	
	param = arg
	
	# Create plots for parameter at once
	layout(matrix(1:4, 2, 2, byrow=T))	
	if (plot.singles == TRUE) {
		for (i in 1:4) {
			if (i == 1) {plot(arg~c(1:length(arg)), main="Trace", type="l")}
			if (i == 2) {densplot(mcmc(arg), main="Density")}
			if (i == 3) {autocorr.plot(arg, auto.layout=F, main="Autocorrelation")}
			if (i == 4) {
				runmean = c()
				for(k in 1:length(arg)) {runmean[k] = mean(arg[1:k])}
				plot(runmean~c(1:length(arg)), main="Running Mean", type="l")
			}
		}
	}	
	
	# Create plots one at a time
	if (plot.singles == TRUE) {
		for (i in 1:4) {
			if (i == 1) {
				plot(arg ~c(1:length(arg)), main="Trace", type="l")
				readline("Press <RETURN> to see plots for next parameter.")
			}
			if (i == 2) {
				densplot(mcmc(arg), main="Density")
				readline("Press <RETURN> to see plots for next parameter.")
			}
			if (i == 3) {
				autocorr.plot(arg, auto.layout=F, main="Autocorrelation")
				readline("Press <RETURN> to see plots for next parameter.")
			}
			if (i == 4) {
				runmean = c()
				for(k in 1:length(arg)) {runmean[k] = mean(arg[1:k])}
				plot(runmean~c(1:length(arg)), main="Running Mean", type="l")
			}
		}
	}	
}

#########################################
### if arg is BayesTraits MCMC object ###
#########################################

if (class(arg) %in% c("MultistateMCMC", "DiscreteMCMC", "ContinuousMCMC", "RegressionMCMC") == T) {	

	# Get parameters to make plots for
	d = names(arg$Results)
	pars = c()
	pars = c(pars, grep("Tree", d))
	pars = c(pars, grep("Lh", d))
	pars = c(pars, grep("q", d))
	pars = c(pars, grep("Alpha", d))
	pars = c(pars, grep("Beta", d))
	if (arg$BasicInfo$Model == "Continuous") {
		if (arg$BasicInfo$Lambda == "ML") {pars = c(pars, grep("Lambda", d))}
		if (arg$BasicInfo$Kappa == "ML") {pars = c(pars, grep("Kappa", d))}
		if (arg$BasicInfo$Delta == "ML") {pars = c(pars, grep("Delta", d))}
		if (arg$BasicInfo$Delta == "ML") {pars = c(pars, grep("OU", d))}
	}
	pars = sort(unique(pars))
	
	# Create plots for each parameter
	layout(matrix(1:4, 2, 2, byrow=T))	
	for (n in 1:length(pars)) {
		plot(arg$Results[,pars[n]]~arg$Results[,1], ylab=d[pars[n]], xlab="Iteration", main=paste("Trace of", d[pars[n]]), type="l")
		densplot(mcmc(arg$Results[,pars[n]]), main=paste("Density of", d[pars[n]]))
		autocorr.plot(arg$Results[,pars[n]], auto.layout=F, main=paste("Autocorrelation", d[pars[n]]))
		runmean = c()
		for(k in 1:length(arg$Results[,pars[n]])) {runmean[k] = mean(arg$Results[,pars[n]][1:k])}
		plot(runmean~arg$Results[,1], ylab=paste("Mean", d[pars[n]]), xlab="Iteration", main=paste("Running Mean of", d[pars[n]]), type="l")
		if (n != length(pars)) {readline("Press <RETURN> to see plots for next parameter.")}
	}
}

}
	




