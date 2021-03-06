###################
### PREPARATION ###
###################

# Make sure 'coda' and 'ape' are installed and loaded
library(coda)
library(ape)

# Navigate to the folder with data and BayesTraits 
setwd('../Data')

# Load functions
for (n in list.files('../R')) {source(paste("../R/", n, sep=""))}

# Load example files
Continuous1 = read.table("Continuous1.txt", header=F, sep="\t")
Continuous2 = read.table("Continuous2.txt", header=F, sep="\t")
Discrete1 = read.table("Discrete1.txt", header=F, sep="\t")
Discrete2 = read.table("Discrete2.txt", header=F, sep="\t")
Tree1 = read.nexus("Tree1.nex")
Tree100 = read.nexus("Tree100.nex")

#####################
### MULTISTATE ML ###
#####################

# MultistateML (all-rates-different)
mml1 = MultistateML(tree=Tree1, data=Discrete1, silent=F)
summary(mml1)

# MultistateML (equal-rates models)
mml2 = MultistateML(tree=Tree1, data=Discrete1, res="q01 q10", silent=F)
summary(mml2)

# Likelihood ratio test (all-rates-different vs. equal-rates models)
lr.test(mml1,mml2)

#######################
### MULTISTATE MCMC ###
#######################

# MultistateMCMC (all-rates-different)
mmc1 = MultistateMCMC(tree=Tree100, data=Discrete1)
summary(mmc1)
plot.btmcmc(mmc1)

# MultistateMCMC (equal-rates models)
mmc2 = MultistateMCMC(tree=Tree100, data=Discrete1, res="q01 q10")
summary(mmc2)
plot.btmcmc(mmc2)

# Bayes factor test (all-rates-different vs. equal-rates models)
bf.test(mmc1,mmc2)

###################
### DISCRETE ML ###
###################

# DiscreteML (independent model)
dml1 = DiscreteML(tree=Tree1, data=Discrete2)
summary(dml1)

# DiscreteML (dependent model)
dml2 = DiscreteML(tree=Tree1, data=Discrete2, dependent=T)
summary(dml2)

# Likelihood-ratio test (independent vs. dependent models)
lr.test(dml1,dml2)

# Visualize rate matrices for ML discrete models
layout(matrix(1:2,1,2))
plot.rates(dml1, main=paste("Indepenent Model: Log likelihood"))
plot.rates(dml2, main=paste("Dependent Model: Log likelihood"))

#####################
### DISCRETE MCMC ###
#####################

# DiscreteMCMC (independent model)
dmc1 = DiscreteMCMC(tree=Tree100, data=Discrete2)
summary(dmc1)
plot.btmcmc(dmc1)

# DiscreteMCMC (dependent model)
dmc2 = DiscreteMCMC(tree=Tree100, data=Discrete2, dependent=T)
summary(dmc2)
plot.btmcmc(dmc2)

# BayesFactor test (indepenent vs. dependent models)
bf.test(dmc1,dmc2)

# Plot rate matrices for MCMC discrete models
layout(matrix(1:2,1,2))
plot.rates(dmc1, main=paste("Indepenent Model: Harmonic Mean"))
plot.rates(dmc2, main=paste("Dependent Model: Harmonic Mean"))

####################################
### CONTINUOUS ML (SINGLE TRAIT) ###
####################################

# ContinuousML (random model)
cml1 = ContinuousML(tree=Tree1, data=Continuous1)
summary(cml1)

# ContinuousML (directional model)
cml2 = ContinuousML(tree=Tree1, data=Continuous1, directional=T)
summary(cml2)

# Likelihood-ratio test (random vs. directional)
lr.test(cml1,cml2)

#####################################
### CONTINUOUS MCMC (CORRELATION) ###
#####################################

# ContinuousMCMC (estimate correlation)
ccmc1 = ContinuousMCMC(tree=Tree100, data=Continuous2, tc = FALSE)
summary(ccmc1)

# ContinuousMCMC (no correlation)
ccmc2 = ContinuousMCMC(tree=Tree100, data=Continuous2, tc = TRUE)
summary(ccmc2)

# BayesFactor test
bf.test(ccmc1,ccmc2)

##########################
### REGRESSION ML/MCMC ###
##########################

# RegressionML 
rml1 = RegressionML(tree=Tree1, data=Continuous2)
summary(rml1)

# RegressionMCMC
rmc1 = RegressionMCMC(tree=Tree100, data=Continuous2)
summary(rmc1)
plot.btmcmc(rmc1)

#########################
### KILL BT PROCESSES ###
#########################

# Start a long MCMC chain, but before MultistateMCMC is finished running
MultistateMCMC(tree=Tree1, data=Discrete1, it=100000000, silent=F)
# press ESC to get the R command prompt back and kill all BayesTraits processes
kill.bt()






