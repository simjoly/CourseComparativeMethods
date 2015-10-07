lr.test = function (model1, model2) {
	Lh1 = model1$Results$Lh
	Lh2 = model2$Results$Lh
	if (class(model1) %in% c("MultistateML", "DiscreteML", "ContinuousML", "RegressionML") == F) {
		stop("LRtest requires two maximum likelihood objects")
	}
	if (class(model2) %in% c("MultistateML", "DiscreteML", "ContinuousML", "RegressionML") == F) {
		stop("LRtest requires two maximum likelihood objects")
	}
	if (length(Lh1) != length(Lh2)) {
		stop("Objects must contain the same number of models.")
	}
	if (length(Lh1[Lh1>Lh2]) != 0 && length(Lh1[Lh1>Lh2]) != length(Lh1)) {
		stop("Models must be nested.")
	}
	if (model1$BasicInfo$Model != model2$BasicInfo$Model) {
		stop("Models must be nested.")
	}
	LRstatistic = c()
	pvalue = c()
	for (n in 1:length(Lh1)) {
		lrs = 2*(Lh1[n] - Lh2[n])
		if (Lh1[n] < Lh2[n]) {lrs = -lrs}
		pv = pchisq(lrs, df=1, lower.tail=F) 
		LRstatistic = c(LRstatistic, lrs)
		pvalue = c(pvalue, pv)
	}
return(data.frame(Lh1, Lh2, LRstatistic, pvalue))
}
