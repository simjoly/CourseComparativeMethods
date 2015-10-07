DiscreteMCMC = function (tree, data, dependent = FALSE, rd = c(), pa = c(), res = c(), resall = c(), mrca = c(), fo = c(), et = c(), it = c(), bi = c(), sa = c(), pr = c(), hp = c(), hpall = c(), rj = c(), rjhp = c(), silent = TRUE) {

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

input = c(3,2)	#Discrete, MCMC
input = c(input, paste("lf ./BTout.log.txt"))	#Logfile
if (length(res) > 0) {for (i in 1:length(res)) {input = c(input, paste("Restrict", res[i]))}}
if (dependent == F) {
	input = c(input, "Restrict q12 q34")
	input = c(input, "Restrict q21 q43")
	input = c(input, "Restrict q13 q24")
	input = c(input, "Restrict q31 q42")
}
if (length(resall) == 1) {input = c(input, paste("RestrictAll", resall))}
if (length(et) > 0) {for (i in 1:length(et)) {input = c(input, paste("et", et[i]))}}
if (length(mrca) > 0) {for (i in 1:length(mrca)) {input = c(input, paste("mrca", paste("mrcaNode", i, sep=""), mrca[i]))}}
if (length(fo) > 0) {for (i in 1:length(fo)) {input = c(input, paste("Fossil", paste("fossilNode", i, sep=""), fo[i]))}}
if (length(it) > 0) {input = c(input, paste("it", format(it, sci=F)))}
if (length(bi) > 0) {input = c(input, paste("bi", format(bi, sci=F)))}
if (length(sa) > 0) {input = c(input, paste("sa", format(sa, sci=F)))}
if (length(rd) > 0) {input = c(input, paste("rd", format(rd, sci=F)))}
if (length(pr) > 0) {for (i in 1:length(pr)) {input = c(input, paste("prior", pr[i]))}}
if (length(pa) == 1) {input = c(input, paste("pa", pa))}
if (length(rj) == 1) {input = c(input, paste("rj", rj))}
if (length(hp) > 0) {for (i in 1:length(hp)) {input = c(input, paste("hp", hp[i]))}}
if (length(hpall) == 1) {input = c(input, paste("Hpall", hpall))}
if (length(rjhp) == 1) {input = c(input, paste("rjhp", rjhp))}
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
Results <- read.table("./BTout.log.txt", skip = (grep("Tree No", FullOut) - 1), sep = "\t", header = TRUE, quote="\"")
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
#system(paste("rm ./BTout.log.txt"))
#system(paste("rm", "./inputfile.txt"))
#system(paste("rm", "./tree.nex"))
#system(paste("rm", "./data.txt"))
#system(paste("rm", "./BTout.log.txt.Schedule.txt"))


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

# Get restrictions
l = TRUE
n = grep("Restrictions", FullOut) + 1
Restrictions = c()
while (l) {
	Restrictions = c(Restrictions, FullOut[n])
	l = (strsplit(FullOut[n+1], split=character(0))[[1]][1] == " ")
	n = n + 1
}
Restrictions = data.frame(Restrictions)

# Create transition matrices
transmatrix = matrix(0, 4, 4, dimnames=list(c("00", "01", "10", "11"), c("00", "01", "10", "11")))
	for (i in 1:4) {
		for (j in 1:4) {
			if (i != j  && sum(i,j) != 5) {
				transmatrix[i,j] = mean(Results[,grep(tail(paste("q", i, j, sep=""), 1), names(Results))])
			} else transmatrix[i,j] = NaN
		}
	}
	
####################
### CREATE BTout ###
####################

# Basic model info
BTout[[1]] = data.frame(Model = "Discrete", Analysis = "MCMC", NumRates = as.numeric(tail(strsplit(FullOut[grep("No of Rates", FullOut)], " ")[[1]], 1)), CharSymbols = tail(strsplit(FullOut[grep("Character Symbols", FullOut)], " ")[[1]], 1), Dependent = dependent, Covarion = FALSE, ReversibleJump = length(which(names(Results)=="Model.string")) > 0)

# MCMC parameters
BTout[[2]] = data.frame(Iterations = as.numeric(tail(strsplit(grep("Iterations", FullOut, value=T), " ")[[1]], 1)), Burnin = as.numeric(tail(strsplit(grep("Burn in", FullOut, value=T), " ")[[1]], 1)), Sample = as.numeric(tail(strsplit(grep("Sample Period", FullOut, value=T), " ")[[1]], 1)), Ratedev = tail(strsplit(grep("Rate Dev", FullOut, value=T), " ")[[1]], 1))

# Priors
BTout[[3]] = PriorIn

# Basic tree info
BTout[[4]] = data.frame(Trees = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+1], " ")[[1]], 1)), Taxa = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+2], " ")[[1]], 1)), Sites = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+3], " ")[[1]], 1)), States = as.numeric(tail(strsplit(FullOut[grep("Tree Information", FullOut)+4], " ")[[1]], 1)))

# Model restrictions
BTout[[5]] = Restrictions

# Fossilized nodes
if (length(fo) > 0) {BTout[[6]] = fo} else BTout[[6]] = "None"

# MRCAs
if (length(mrca) > 0) {BTout[[7]] = mrca} else BTout[[7]] = "None"

# Excluded taxa
if (length(et) > 0) {BTout[[8]] = et} else BTout[[8]] = "None"

# Rate matrices
BTout[[9]] = list(transmatrix)

# Results summary
BTout[[10]] = Results

# Phylogeny used
BTout[[11]] = tree

# Data used
BTout[[12]] = data.frame(Data = data)

# Schedule
BTout[[13]] = Schedule

# Add names to object BTout
names(BTout) = c("BasicInfo", "MCMC", "PriorInfo", "TreeInfo", "Restrictions", "Fossil", "AddMRCA", "ExcludedTaxa", "RateMatrix", "Results", "Tree", "Data", Rate)

# Add names to object BTout
class(BTout) = "DiscreteMCMC"
return(BTout)

}


