kill.bt = function() {
	if (Sys.info()[['sysname']]=="Windows") { 		# Windows
		if (Sys.getenv("R_ARCH")=="/i386") { 		# Running 32 bit version of windows
			shell("Taskkill /fi \"Imagename eq BayesTraitsV2beta_win32.exe\"")
		} else if (Sys.getenv("R_ARCH")=="/x64") { 	# Running 64 bit version of windows
			shell("Taskkill /fi \"Imagename eq BayesTraitsV2beta_win64.exe\"")
		} else {
			stop("Cannot determine Windows architecture (32 or 64 bits)")
		}
	} else if (Sys.info()[['sysname']]=="Darwin") { # MAC
		jobs = system("pgrep BayesTraitsV2_mac", intern=T)
		for (n in 1:length(jobs)) {system(paste("kill", jobs[n]))}
	} else { 										# Linux
		jobs = system("pgrep BayesTraitsV2_linux", intern=T)
		for (n in 1:length(jobs)) {system(paste("kill", jobs[n]))}
	}
}
