# Last update by Randi Griffin on June 4, 2014

MultistateML = function (tree, data, mlt = c(), res = c(), resall = c(), mrca = c(), fo = c(), et = c(), silent = TRUE) {
	
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

#########################################
### CREATE TREE, DATA and INPUT FILES ###
#########################################

input = c(1,1)	
input = c(input, paste("lf ./BTout.log.txt"))
if (length(mlt) > 0) {input = c(input, paste("mlt", as.numeric(mlt)))}
if (length(res) > 0) {for (i in 1:length(res)) {input = c(input, paste("Restrict", res[i]))}}
if (length(resall) == 1) {input = c(input, paste("RestrictAll", resall))}
if (length(mrca) > 0) {for (i in 1:length(mrca)) {input = c(input, paste("mrca", paste("mrcaNode", i, sep=""), mrca[i]))}}
if (length(fo) > 0) {for (i in 1:length(fo)) {input = c(input, paste("Fossil", paste("fossilNode", i, sep=""), fo[i]))}}
if (length(et) > 0) {for (i in 1:length(et)) {input = c(input, paste("et", et[i]))}}
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

# Get states
options(warn=-1)
rcs = tail(strsplit(FullOut[grep("Character Symbols", FullOut)], " ")[[1]],1)
states = substr(rcs, seq(1,nchar(rcs),2), seq(1,nchar(rcs),2))
options(warn=0)	

# Create transition matrices
Matrix = list()
for (n in 1:nrow(Results)) {
	mat = matrix(0, length(states), length(states), dimnames=list(as.character(states), as.character(states)))
	for (i in 1:length(states)) {
		for (j in 1:length(states)) {
			if (i != j) {
				mat[i,j] = Results[n,grep(tail(paste("q", states[i], states[j], sep=""), 1), names(Results))]
			} else mat[i,j] = NaN
		}
	}
	Matrix[[n]] = mat
}

####################
### CREATE BTout ###
####################

# Basic model info
BTout[[1]] = data.frame(Model = "Multistates", Analysis = "Maximum likelihood", MLattempts = as.numeric(tail(strsplit(FullOut[grep("ML attempt", FullOut)], " ")[[1]], 1)), NumRates = as.numeric(tail(strsplit(FullOut[grep("No of Rates", FullOut)], " ")[[1]], 1)), CharSymbols = tail(strsplit(FullOut[grep("Character Symbols", FullOut)], " ")[[1]], 1), Covarion = FALSE)

# Basic tree info
BTout[[2]] = data.frame(Trees = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+1], " ")[[1]], 1)), Taxa = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+2], " ")[[1]], 1)), Sites = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+3], " ")[[1]], 1)), States = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+4], " ")[[1]], 1)))

# Model restrictions
BTout[[3]] = Restrictions

# Fossilized nodes
if (length(fo) > 0) {BTout[[4]] = fo} else BTout[[4]] = "None"

# MRCAs
if (length(mrca) > 0) {BTout[[5]] = mrca} else BTout[[5]] = "None"

# Excluded taxa
if (length(et) > 0) {BTout[[6]] = et} else BTout[[6]] = "None"

# Rate matrices
BTout[[7]] = Matrix

# Results summary
BTout[[8]] = Results

# AIC
BTout[[9]] = c(AIC = ((-2)*BTout[[8]]$Lh)+(2*NumRateEst), NumPar = NumRateEst)

# Phylogeny used
BTout[[10]] = tree

# Data used
BTout[[11]] = data.frame(Data = data)

# Add names to objects in BTout
names(BTout) = c("BasicInfo", "TreeInfo", "Restrictions", "Fossil", "AddMRCA", "ExcludedTaxa", "RateMatrix", "Results", "AIC", "Tree", "Data")

# Return object of class MultistateML
class(BTout) = "MultistateML"
return(BTout)

}




