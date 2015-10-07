ContinuousMCMC = function(tree, data, directional = FALSE, lambda = 1, kappa = 1, delta = 1, ou = 0, tc = FALSE, it = c(), bi = c(), sa = c(), rd = c(), silent = TRUE) {

###################################
### CHECK FOR REQUIRED PACKAGES ###
###################################

if (("ape" %in% rownames(installed.packages())) == F) {stop("You must install the package 'ape'.")}
library(ape)

#############################################
### CHECK FOR MISMATCHES IN SPECIES NAMES ###
#############################################

species = c()
if(class(tree)=="phylo") {species=tree$tip.label} 
if(class(tree)=="multiPhylo") {species=tree$tip.label$tip.label}
MissingFromData = setdiff(species, data[,1])
MissingFromTree = setdiff(data[,1], species)

if (length(MissingFromData)>0) {
	warning(paste("No match found in the data:", paste(MissingFromData, collapse=", ")))
}
if (length(MissingFromTree)>0) {
	warning(paste("No match found in the phylogeny:", paste(MissingFromTree, collapse=", ")))
}
if(length(MissingFromData)>0 | length(MissingFromTree)>0) {
	stop("Species in your phylogeny and data must match up exactly.")
}

#########################
### CREATE INPUT FILE ###
#########################

if (directional == F) {dir = 4} else {dir = 5}
if (lambda == "ML") {L = ""} else {L = lambda}
if (kappa == "ML") {K = ""} else {K = kappa}
if (delta == "ML") {D = ""} else {D = delta}
if (ou == "ML") {O = ""} else {O = ou}
input = c(dir, 2)
input = c(input, paste("lf ./BTout.log.txt"))	
input = c(input, paste("lambda", L))
input = c(input, paste("kappa", K))
input = c(input, paste("delta", D))
input = c(input, paste("ou", O))
if (tc == TRUE) {input = c(input, "tc")}
if (length(it) > 0) {input = c(input, paste("it", format(it, sci=F)))}
if (length(bi) > 0) {input = c(input, paste("bi", format(bi, sci=F)))}
if (length(sa) > 0) {input = c(input, paste("sa", format(sa, sci=F)))}
if (length(rd) > 0) {input = c(input, paste("rd", format(rd, sci=F)))}
input = c(input, "Schedule True")
input = c(input, "run")	
write(input, file="./inputfile.txt") 
write.nexus(tree, file="./tree.nex", translate=T)	
write.table(data, file="./data.txt", quote=F, col.names=F, row.names=F)

#########################
### RUN & GET RESULTS ###
#########################
	
# Run
if (Sys.info()[['sysname']]=="Windows") { 		# Windows
	if (Sys.getenv("R_ARCH")=="/i386") { 		# Running 32 bit version of windows
		shell(paste("BayesTraitsV2beta_win32.exe", "tree.nex", "data.txt", "< inputfile.txt"), ignore.stdout = silent)
	} else if (Sys.getenv("R_ARCH")=="/x64") { 	# Running 64 bit version of windows
		shell(paste("BayesTraitsV2_win64.exe", "tree.nex", "data.txt", "< inputfile.txt"), ignore.stdout = silent)
	} else {
		stop("Cannot determine Windows architecture (32 or 64 bits)")
	}
} else if (Sys.info()[['sysname']]=="Darwin") { # MAC
   system(paste("./BayesTraitsV2_mac", "./tree.nex", "./data.txt", "< ./inputfile.txt"), ignore.stdout = silent)
} else { 										# Linux
   system(paste("./BayesTraitsV2_linux", "./tree.nex", "./data.txt", "< ./inputfile.txt"), ignore.stdout = silent)
}

# Get BT output
FullOut = scan(file = "./BTout.log.txt", what="c", quiet=T, sep="\n")
Results <- read.table("./BTout.log.txt", skip = (grep("Tree No", FullOut) - 1), sep = "\t", header = TRUE)
Results = Results[,-ncol(Results)]	

# Get BT schedule
RateTreeMove = scan(file="./BTout.log.txt.Schedule.txt", what="c", nlines=2, sep="\t", quiet=T)
Skip = length(RateTreeMove)/2
Schedule = read.table("./BTout.log.txt.Schedule.txt", skip= Skip, header=T, sep="\t")
Schedule = Schedule[,-ncol(Schedule)]
Rate = paste(RateTreeMove, collapse=" ")

# Delete files from disk
if (Sys.info()[['sysname']]=="Windows") { #Windows
	shell(paste("del BTout.log.txt inputfile.txt tree.nex data.txt BTout.log.txt.Schedule.txt"))
}	else {
	system(paste("rm ./BTout.log.txt ./inputfile.txt ./tree.nex ./data.txt ./BTout.log.txt.Schedule.txt"))
}

######################
### PROCESS OUTPUT ###
######################

BTout = list()

# Get prior information
l = TRUE
n = grep("Prior Information", FullOut) + 1
PriorIn = c()
while (l) {
	if ((strsplit(FullOut[n+1], split=character(0))[[1]][1] == " ")) {
		PriorIn = c(PriorIn, FullOut[n])
		n = n + 1
		} else l = FALSE
}
PriorIn = data.frame(PriorIn)

# Get number of estimated rates and restrictions
l = TRUE
n = grep("Restrictions", FullOut) + 1
NumRateEst = 0
Restrictions = c()
while (l) {
	Restrictions = c(Restrictions, FullOut[n])
	if (tail(strsplit(FullOut[n], " ")[[1]], 1) == "None") {NumRateEst = NumRateEst + 1}
	l = (strsplit(FullOut[n+1], split=character(0))[[1]][1] == " ")
	n = n + 1
}
Restrictions = data.frame(Restrictions)

####################
### CREATE BTout ###
####################

# Basic model info
BTout[[1]] = data.frame(Model = "Continuous", Analysis = "MCMC", NumRates = as.numeric(tail(strsplit(FullOut[grep("No of Rates", FullOut)], " ")[[1]], 1)), Lambda = lambda, Kappa = kappa, Delta = delta, OU = ou, Directional = directional, TestCorrel = tc)

# MCMC parameters
BTout[[2]] = data.frame(Iterations = as.numeric(tail(strsplit(grep("Iterations", FullOut, value=T), " ")[[1]], 1)), Burnin = as.numeric(tail(strsplit(grep("Burn in", FullOut, value=T), " ")[[1]], 1)), Sample = as.numeric(tail(strsplit(grep("Sample Period", FullOut, value=T), " ")[[1]], 1)), Ratedev = tail(strsplit(grep("Rate Dev", FullOut, value=T), " ")[[1]], 1))

# Priors
BTout[[3]] = PriorIn

# Basic tree info
BTout[[4]] = data.frame(Trees = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+1], " ")[[1]], 1)), Taxa = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+2], " ")[[1]], 1)), Sites = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+3], " ")[[1]], 1)))

# Model restrictions
BTout[[5]] = Restrictions

# Results summary
BTout[[6]] = Results

# Phylogeny used
BTout[[7]] = tree

# Data used
BTout[[8]] = data.frame(Data = data)

# Schedule
BTout[[9]] = Schedule

# Add names to object BTout
names(BTout) = c("BasicInfo", "MCMC", "PriorInfo", "TreeInfo", "Restrictions", "Results", "Tree", "Data", Rate)

# Return object of class ContinuousMCMC
class(BTout) = "ContinuousMCMC"
return(BTout)

}


