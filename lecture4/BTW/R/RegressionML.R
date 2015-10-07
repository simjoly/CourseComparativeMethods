RegressionML = function(tree, data, lambda = 1, kappa = 1, delta = 1, ou = 1, mlt = c(), silent=TRUE) {

###################################
### CHECK FOR REQUIRED PACKAGES ###
###################################

if (("ape" %in% rownames(installed.packages())) == F) {stop("You must install the package 'ape'.")}
library(ape)

########################################
### REMOVE SPACES FROM SPECIES NAMES ###
########################################

tree$tip.label = gsub(" ", "_", tree$tip.label)
data[,1] = gsub(" ", "_", data[,1])

#############################################
### CHECK FOR MISMATCHES IN SPECIES NAMES ###
#############################################

MissingFromData = setdiff(tree$tip.label, data[,1])
MissingFromTree = setdiff(data[,1], tree$tip.label)

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

input = c(6, 1)
if (lambda == "ML") {L = ""} else {L = lambda}
if (kappa == "ML") {K = ""} else {K = kappa}
if (delta == "ML") {D = ""} else {D = delta}
if (ou == "ML") {O = ""} else {O = ou}
input = c(input, paste("lf ./BTout.log.txt"))	
input = c(input, paste("lambda", L))
input = c(input, paste("kappa", K))
input = c(input, paste("delta", D))
input = c(input, paste("ou", O))
if (length(mlt) > 0) {input = c(input, paste("mlt", as.numeric(mlt)))}
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
Results = read.table("./BTout.log.txt", skip = (grep("Tree No", FullOut) - 1), sep = "\t", header = TRUE)
Results = Results[,-ncol(Results)]

# Delete files from disk
if (Sys.info()[['sysname']]=="Windows") { #Windows
	shell(paste("del BTout.log.txt inputfile.txt tree.nex data.txt"))
}	else {
	system(paste("rm ./BTout.log.txt ./inputfile.txt ./tree.nex ./data.txt"))
}

######################
### PROCESS OUTPUT ###
######################

BTout = list()

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
BTout[[1]] = data.frame(Model = "Continuous", Analysis = "Maximum likelihood", MLattempts = as.numeric(tail(strsplit(FullOut[grep("ML attempt", FullOut)], " ")[[1]], 1)), NumRates = as.numeric(tail(strsplit(FullOut[grep("No of Rates", FullOut)], " ")[[1]], 1)), Lambda = lambda, Kappa = kappa, Delta = delta, OU = ou)

# Basic tree info
BTout[[2]] = data.frame(Trees = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+1], " ")[[1]], 1)), Taxa = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+2], " ")[[1]], 1)), Sites = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+3], " ")[[1]], 1)))

# Model restrictions
BTout[[3]] = Restrictions

# Results summary
BTout[[4]] = Results

# AIC
BTout[[5]] = c(AIC = ((-2)*BTout[[4]]$Lh)+(2*BTout[[1]]$NumRates), NumPar = BTout[[1]]$NumRates)

# Phylogeny used
BTout[[6]] = tree

# Data used
BTout[[7]] = data.frame(Data = data)

# Add names to objects in BTout
names(BTout) = c("BasicInfo", "TreeInfo", "Restrictions", "Results", "AIC", "Tree", "Data")

# Return object of class RegressionML
class(BTout) = "RegressionML"
return(BTout)

}

