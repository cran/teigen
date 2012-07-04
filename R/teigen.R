teigen <-
function(x, Gs=1:9, models="all", init="kmeans", scale=TRUE, dfstart=50, clas=0, known=NULL, training=NULL, gauss=FALSE, dfupdate=TRUE, eps=0.1){
	origx <- x
	teigenModels <- list()
	teigenModels[["multivariate"]] <- c("UUUU", "UUUC","CUCU", "CUCC","CUUU", "CUUC","CCCU","CCCC",
										"CIUU", "CIUC","CICU", "CICC","UIIU","UIIC","CIIU","CIIC",
										"UIUU","UIUC","UCCU", "UCCC","UUCU", "UUCC","UICU", "UICC")
	teigenModels[["univariate"]] <- c("univUU", "univUC", "univCU", "univCC")
	teigenModels[["dfconstrained"]] <- c("UUUC","CUCC","CUUC","CCCC",
										 "CIUC", "CICC","UIIC","CIIC",
										 "UIUC","UCCC","UUCC", "UICC")
	teigenModels[["dfunconstrained"]] <- c("UUUU","CUCU","CUUU","CCCU",
										   "CIUU", "CICU","UIIU","CIIU",
										   "UIUU","UCCU","UUCU", "UICU")
	noconv <- FALSE
	if(is.vector(x)){
		x <- matrix(origx,length(origx),1)
		univar <- TRUE
	}
	else{
		if(ncol(x)<2){
			univar <- TRUE
		}
		else{
			if(length(ncol(x))<0){
				univar <- TRUE
			}
			else{
				univar <-  FALSE
			}
		}
	}
	if(univar){
		if(length(models)>1){
			if(gauss){
				dfstart <- 200
				gauss <- TRUE
				dfupdate <- FALSE
				models <- c("univUC","univCC")
			}
			else{
				if(any(!models%in%teigenModels$univariate)){
					models <- teigenModels$univariate
				}
			}
		}
		else{
			if(models=="mclust"|models=="gaussian"){
				dfstart <- 200
				gauss <- TRUE
				dfupdate <- FALSE
				models <- c("univUC","univCC")
			}
			else{
				if((!models%in%teigenModels$univariate) & gauss){
					models <- teigenModels$univariate
				}
			}
		}
	}
	known <- factor(known)
	if(nrow(x)<ncol(x)){
		warning("Dimensionality of data exceeds the number of samples, spherical models suggested (those with 'I' as a second letter)") 
	}
	for(i in 1:length(models)){
		if(!any(models[i]==c("UUUU", "UUUC","CUCU", "CUCC","CUUU", "CUUC","CCCU","CCCC",
					 "CIUU", "CIUC","CICU", "CICC","UIIU","UIIC","CIIU","CIIC",
					 "UIUU","UIUC","UCCU", "UCCC","UUCU", "UUCC","UICU", "UICC",
					 "univUU", "univUC", "univCU", "univCC","all","gaussian","mclust",
					 "dfconstrained","dfunconstrained","univariate"))){
			stop("You have specified at least one unknown model abbreviation...please select a different set of models")
			return(NULL)
		}
	}
	if(gauss){
		dfstart <- 200
		dfupdate <- FALSE
	}
	if(scale){
		x <- scale(x, center=TRUE, scale=TRUE)
	}
	p <- ncol(x)
	n <- nrow(x)
	if(length(models)==1){
		if(models=="dfunconstrained"){
			models <- teigenModels$dfunconstrained
		}
		else{
			if(models=="all"){
				if(univar){
					models <- teigenModels$univariate
				}
				else{
					models <- teigenModels$multivariate
				}
			}
			else{
				if(models=="gaussian"){
					models <- teigenModels$dfconstrained
					dfstart <- 200
					gauss <- TRUE
					dfupdate <- FALSE
				}
				else{
					if(models=="mclust"){
						models <- c("UUUC","CUCC","CCCC","CIUC","CICC","UIIC","CIIC","UIUC","UUCC","UICC")
						dfstart <- 200
						gauss <- TRUE
						dfupdate <- FALSE
					}
					else{
						if(models=="dfconstrained"){
							models <- teigenModels$dfconstrained
						}
						else{
							if(models=="univariate"){
								models <- teigenModels$univariate
							}
						}
					}
				}
			}
		}
	}
	if(any(c("UCUU","UCUC","CCUU","CCUC") %in% models)){
		stop("Models UCCU, UCCC, CCUU ,CCUC are not yet available...please select a different set of models")
		return(NULL)
	}
	hh8 <- dfupdate
	zlist3 <- list()
	dff <- list()
	it <- list()
	store <- list()
	meanlist <- list()
	lambdalist <- list()
	dlist <- list()
	alist <- list()
	siglist <- list()
	if(clas>0){
		if(length(known)!=n){
			stop("Known classifications vector not given, or not the same length as the number of samples (see help file)")
			return(NULL)
		}
		testindex <- sample(1:n, ceiling(n*(clas/100)))
		kno <- vector(mode="numeric", length=n)
		kno[testindex] <- 1
		unkno <- (kno-1)*(-1)
		Gs <- length(unique(known))
	}
	if(length(training)>0){
		if(length(known)!=n){
			stop("Known classifications vector not given, or not the same length as the number of samples (see help file)")
			return(NULL)
		}
		testindex <- training
		kno <- vector(mode="numeric", length=n)
		kno[testindex] <- 1
		unkno <- (kno-1)*(-1)
		Gs <- length(unique(known[training]))
		clas <- length(training)/nrow(x)
	}
	gvec <- 1:max(Gs)
	gstuff <- paste("G=",gvec,sep="")
	bic <- matrix(-Inf, length(models), max(Gs))
	icl <- matrix(-Inf, length(models), max(Gs))
	logls <- matrix(-Inf, length(models), max(Gs))
	unc <- matrix(Inf,length(models),max(Gs))
	zmatin <- list()
	for(G in Gs){
		if(G==1){
			zmatin[[G]] <- matrix(1,n,1)
		}
		else{
			if(is.character(init)){
				if(init == "disc" | init == "hard"){
					zmatin[[G]] <- tdiscrandz(n,G)
				}
				if(init == "cont" | init == "soft"){
					zmatin[[G]] <- tcontrandz(n,G)
				}
				if(init == "uniform"){
					if(clas>0){
						zmatin[[G]] <- tuniformz(n,G,clas,kno,known)
					}
					else{
						stop("Uniform initialization not available for clustering.")
						return(NULL)
					}
				}
				if(init == "kmeans"){
					zmatin[[G]] <- tkmeansz(x,n,G,known,kno,testindex,clas)
				}
			}
			else{
				zmatin[[G]] <- tgivenz(n,G,known,init[[G]],testindex,clas)
			}
		}
	}
for(modnum in 1:length(models)){
	mod <- models[modnum]
	submod13 <- substring(mod,1,3) 
	zlist2 <- list()
	it2 <- list()
	df2 <- list()
	siglist2 <- list()
	meanlist2 <- list()
	lambdalist2 <- list()
	dlist2 <- list()
	alist2 <- list()
	for(G in Gs){ 
	singular <- 0
		killit <- FALSE
		if(G==1){
			if(length(models)>1){
				CCCCgroup <- c("CCCC", "CCCU", "CUCC", "CUCU", "CUUC","CUUU","UCCC","UCCU","UUCU","UUCC", "UUUC","UUUU")
				if(any(mod==CCCCgroup)){
					cccdum <- models[models %in% CCCCgroup]
					if(length(cccdum)>0){
						if(mod!=cccdum[1]){
							killit <- TRUE
						}
					}
				}
				CICCgroup <- c("CICC","CICU","UICC","UICU","CIUC","CIUU","UIUC","UIUU")
				if(any(mod==CICCgroup)){
					cicdum <- models[models %in% CICCgroup]
					if(length(cicdum)>0){
						if(mod!=cicdum[1]){
							killit <- TRUE
						}
					}
				}
				CIICgroup <- c("CIIC", "CIIU", "UIIC","UIIU")
				if(any(mod==CIICgroup)){
					ciidum <- models[models %in% CIICgroup]
					if(length(ciidum)>0){
						if(mod!=ciidum[1]){
							killit <- TRUE
						}
					}
				}
				univgroup <- c("univCC","univCU","univUC","univUU")
				if(any(mod==univgroup)){
					unidum <- models[models %in% univgroup]
					if(length(unidum)>0){
						if(mod!=unidum[1]){
							killit <- TRUE
						}
					}
				}
			}
		}
	zmat <- zmatin[[G]]
	vg <- tvginit(dfstart,G)
	ng <- tngupdate(zmat)
	if(any(ng<1.5)){break}
	pig <- tpigupdate(ng,n)
		mug <- matrix(0,G,p)
		sg <- array(0,dim=c(p,p,G))
		for(g in 1:G){
			wtcov <- cov.wt(x,wt=zmat[,g],method="ML")
			mug[g,] <- wtcov$center
			sg[,,g] <- wtcov$cov
		}
	sgc <- tsginitc(G,sg,pig,p,n,x)
	if(any(submod13==c("CCC","CIC")) | mod=="univCC" | mod=="univCU"){
		sg[,,] <- sgc 
	}
	if(!univar){
			for(g in 1:G){
				test <- rcond(sg[,,g])
					if(test <= sqrt(.Machine$double.eps)){
							singular <- 1
					}
			}
			if(singular==1){
				break
			}
		if(all(submod13!=c("CCC","UUU"))){
			if(any(submod13==c("UCC","CCU","UCU","UUC","UIC","UCU"))){
				lambdag <- tlambdaginit(p,G,sg,mod)
				dg <- tdginit(p,G,sg,mod,sgc)
				ag <- taginit(p,G,sg,mod,sgc)
			}
			else{
				lambdag <- tlambdagupdate(G,mod,sg,dg,ag,p,n,ng,submod13)
				dg <- tdgupdate(p,G,mod,sg,lambdag,ng,ag,submod13)
				if(any(is.logical(dg))){
					killit <- TRUE
				}
				else{
					ag <- tagupdate(p,G,mod,sg,lambdag,ng,n,dg,submod13)
					if(any(is.logical(ag))){
						break
					}
				}
			}
		}
		else{
			dg <- Inf
			ag <- Inf
			lambdag <- Inf
		}
	}
	if(!killit){
		sigma <- tsigmaup(p,G,sg,lambdag,dg,ag,mod,univar,submod13)
		if(!univar){
				for(g in 1:G){
					test <- rcond(sigma[,,g])
						if(test <= sqrt(.Machine$double.eps)){
								singular <- 1
						}
				}
				if(singular==1){break}
			sigmainv <- tsigmainvup(p,G,sigma)
		}
		w <- twinit(x,n,G,mug,sigmainv,vg,p,sg,zmat,univar,sigma)
	}
		cycle <- 0
		dhfgs78 <- vg
		conv <- 0
		num <- matrix(0,n,G)
		ft <- matrix(0,n,G)
		logl <- NaN
	while(conv != 1){
		if(killit){break}
			ng <- tngupdate(zmat)
			if(any(ng<1.5)){break}
			pig <- tpigupdate(ng,n)
			mug <- tmugupdate(G,zmat,w,x,p,univar)
			if(hh8){
				testing <- try(jk861 <- yxf7(mod,dhfgs78,ng,zmat,w,G,p,n,x,mug,sigmainv),silent=TRUE)
				if(all(is.finite(testing))){	
					dhfgs78 <- jk861
				}
				else{break}
			}
				delta <- deltaup(x,mug,sigma,sigmainv,G,n,univar)
				zmat <- tzupdate(x,G,pig,dhfgs78,p,mug,sigmainv,n,sigma,clas,kno,known,unkno,univar,delta)
				if(any(is.nan(zmat))){
					break
				}
				w <- twupdate(x,n,G,mug,sigmainv,dhfgs78,p,univar,sigma,delta)
			ng <- tngupdate(zmat)
			if(any(ng<1.5)){break}
			sg <- tsgupdate(p,G,n,x,mug,zmat,w,ng,mod,pig,submod13)
		if(!univar){
			if(all(submod13!=c("CCC","UUU"))){
				if(any(submod13==c("UCC","UUC","UIC"))){
					mcyc <- 0
					fmin <- 0
					mtest <- Inf
					while(mtest > 0.05){
						mcyc <- mcyc + 1
						lambdag <- tlambdagupdate(G,mod,sg,dg,ag,p,n,ng,submod13)
						dg <- tdgupdate(p,G,mod,sg,lambdag,ng,ag,submod13)
						if(any(is.logical(dg))){
							killit <- TRUE
							break
						}
						ag <- tagupdate(p,G,mod,sg,lambdag,ng,n,dg,submod13)
						if(any(is.logical(ag))){
							killit <- TRUE
							break
						}
						fmin[mcyc] <- tfminup(mod,G,sg,dg,ag,p,ng,lambdag,submod13)
						if(mcyc > 1){
							mtest <- fmin[mcyc-1] - fmin[mcyc] 
						}
					}
				}
				else{
						lambdag <- tlambdagupdate(G,mod,sg,dg,ag,p,n,ng,submod13)
						dg <- tdgupdate(p,G,mod,sg,lambdag,ng,ag,submod13)
					if(any(is.logical(dg))){
						break
					}
						ag <- tagupdate(p,G,mod,sg,lambdag,ng,n,dg,submod13)
					if(any(is.logical(ag))){
						break
					}
				}
				if(killit){break}
			}
		}
			sigma <- tsigmaup(p,G,sg,lambdag,dg,ag,mod,univar,submod13)
		if(!univar){
				for(g in 1:G){
					test <- rcond(sigma[,,g])
					if(test <= sqrt(.Machine$double.eps)){
							singular <- 1
					}
				}
				if(singular==1){break}
			sigmainv <- tsigmainvup(p,G,sigma)
		}
			delta <- deltaup(x,mug,sigma,sigmainv,G,n,univar)
			ft <- exp(tft(x,G,pig,dhfgs78,p,mug,sigmainv,n,sigma,clas,kno,known,unkno,univar,delta))
			zmat <- ft/rowSums(ft)
		if(clas>0){
			zmat <- unkno*zmat
			for(i in 1:n){
				if(kno[i]==1){
					zmat[i, known[i]] <- 1
				}
			}
		}
		if(any(is.nan(zmat))){
			break
		}
			w <- twupdate(x,n,G,mug,sigmainv,dhfgs78,p,univar,sigma,delta)
			ng <- tngupdate(zmat)
			if(any(ng<1.5)){break}
			cycle <- cycle + 1
				logl[cycle]<- sum(log(rowSums(ft)))
			if(is.na(logl[cycle])){break}
		if(cycle>3){
				ak <- (logl[cycle]-logl[cycle-1])/(logl[cycle-1]-logl[cycle-2])
				linf <- logl[cycle-1] + (logl[cycle]-logl[cycle-1])/(1-ak)
				if(abs(linf-logl[cycle-1]) < eps){
					conv<-1
				}
				if((logl[cycle]-logl[cycle-1])<0){
					break
				}
		}
}
	if(conv==1){
		logls[modnum,G] <- max(logl)
		bic[modnum,G] <- tBICcalc(conv,G,p,mod,logl,n,gauss,univar,submod13)
		icl[modnum,G] <- tICLcalc(conv,n,zmat,bic,modnum,G)
		it2[[G]] <- cycle
		df2[[G]] <- dhfgs78
		meanlist2[[G]] <- mug
		if(!univar){
			lambdalist2[[G]] <- lambdag
			dlist2[[G]] <- dg
			alist2[[G]]<- ag
		}
		siglist2[[G]] <- sigma
		zlist2[[G]] <- zmat
	}else{
		if(G>1){
				noconv <- TRUE
		}
	}
		}
		it[[modnum]] <- it2
		zlist3[[modnum]] <- zlist2
		dff[[modnum]] <- df2
		meanlist[[modnum]] <- meanlist2
		if(!univar){
			lambdalist[[modnum]] <- lambdalist2
			dlist[[modnum]] <- dlist2
			alist[[modnum]] <- alist2
		}
		siglist[[modnum]] <- siglist2
	}
	if(all(is.infinite(bic))){
		stop("No models converged. Try different models, number of components, or different initialization.")
		return(NULL)
	}
	dimnames(bic) <- list(models,gstuff[1:max(Gs)])
	dimnames(logls) <- list(models,gstuff[1:max(Gs)])
	dimnames(icl) <- list(models,gstuff[1:max(Gs)])
	dimnames(unc) <- list(models,gstuff[1:max(Gs)])
	unc[,1] <- Inf
	maxes <- which(bic==max(bic), arr.ind=TRUE)
	maxicl <- which(icl==max(icl), arr.ind=TRUE)
	minunc <- which(unc==min(unc), arr.ind=TRUE)
	known <- as.character(known)
	known[is.na(known)] <- "unknown"
	known <- factor(as.character(known))
	if(nrow(maxes)>1){
		warning("Maximum BIC tie between two or more models")
		bestmodnum <- maxes[1:nrow(maxes),1]
		bestmod <- models[bestmodnum]
		bestg <- maxes[1:nrow(maxes),2]
		itf <- "MULTIPLE"
		dff1 <- "MULTIPLE" 
		bestz <- "MULTIPLE"
		bestzmap <- "MULTIPLE"
		adjrand <- "MULTIPLE"
		tab <- "MULTIPLE"
		bestmean <- "MULTIPLE"
		bestlambda <- "MULTIPLE"
		bestd <- "MULTIPLE"
		besta <- "MULTIPLE"
		bestsig <- "MULTIPLE"
	}
	if(nrow(maxes)==1){
		bestmodnum <- maxes[1]
		bestmod <- models[bestmodnum]
		bestg <- maxes[2]
		bestz <- zlist3[[bestmodnum]][[bestg]]
		dff1 <- dff[[bestmodnum]][[bestg]]
		itf <- it[[bestmodnum]][[bestg]]
		bestzmap <- apply(bestz,1,which.max)
		if(clas>0){
			newmap <- bestzmap
			newmap[testindex] <- NA
			newknown <- known
			newknown[testindex] <- NA
			tab <- table(known,newmap)
		}
		else{
			if(length(known)>0){
				tab <- table(known,bestzmap)
			}
		}
		bestmean <- meanlist[[bestmodnum]][[bestg]]
		if(!univar){
			if(models[bestmodnum]%in%c("UUUU","UUUC","CCCC","CCCU")){
				decom <- list()
				bestd <- array(0, dim=c(p, p, bestg))
				bestlambda <- NA
				besta <- bestd
				for(g in 1:bestg){
					decom[[g]] <- eigen(siglist[[bestmodnum]][[bestg]][,,g])
					bestd[,,g] <- decom[[g]]$vectors
					eigvals <- decom[[g]]$values
					bestlambda[g] <- prod(eigvals)^(1/p)
					besta[,,g] <- diag(eigvals)/bestlambda[g]
				}
			}
			else{
				bestlambda <- lambdalist[[bestmodnum]][[bestg]]
				bestd <- dlist[[bestmodnum]][[bestg]]
				besta <- alist[[bestmodnum]][[bestg]]
			}
		}
		bestsig <- siglist[[bestmodnum]][[bestg]]
	}
	if(nrow(maxicl)>1){
		warning("Maximum ICL tie between two or more models")
		bestmodnumicl <- maxicl[1:nrow(maxicl),1]
		bestmodicl <- models[bestmodnumicl]
		bestgicl <- maxicl[1:nrow(maxicl),2]
		dff1icl <- "MULTIPLE"
		bestzicl <- "MULTIPLE"
		bestzmapicl <- "MULTIPLE"
		adjrandicl <- "MULTIPLE"
		itficl <- "MULTIPLE"
		tabicl <- "MULTIPLE"
		bestmeanicl <- "MULTIPLE"
		bestlambdaicl <- "MULTIPLE"
		bestdicl <- "MULTIPLE"
		bestaicl <- "MULTIPLE"
		bestsigicl <- "MULTIPLE"
	}
	if(nrow(maxicl)==1){
		bestmodnumicl <- maxicl[1]
		bestmodicl <- models[bestmodnumicl]
		bestgicl <- maxicl[2]
		bestzicl <- zlist3[[bestmodnumicl]][[bestgicl]]
		dff1icl <- dff[[bestmodnumicl]][[bestgicl]]
		itficl <- it[[bestmodnumicl]][[bestgicl]]
		bestzmapicl <- apply(bestzicl,1,which.max)
		if(clas>0){
			newmapicl <- bestzmapicl
			newmapicl[testindex] <- NA
			newknown <- known
			newknown[testindex] <- NA
			tabicl <- table(known,newmapicl)
		}
		else{
			if(length(known)>0){
				tabicl <- table(known,bestzmapicl)
			}
		}
		bestmeanicl <- meanlist[[bestmodnumicl]][[bestgicl]]
		if(!univar){
			bestlambdaicl <- lambdalist[[bestmodnumicl]][[bestgicl]]
			bestdicl <- dlist[[bestmodnumicl]][[bestgicl]]
			bestaicl <- alist[[bestmodnumicl]][[bestgicl]]
		}
		bestsigicl <- siglist[[bestmodnumicl]][[bestgicl]]
	}
	icllist <- list()
	parameters <- list()
	parametersicl <- list()
	parameters[["mean"]] <- bestmean
	if(!univar){
		parameters[["a"]] <- besta
		parameters[["d"]] <- bestd
		parameters[["lambda"]] <- bestlambda
		parametersicl[["a"]] <- bestaicl
		parametersicl[["d"]] <- bestdicl
		parametersicl[["lambda"]] <- bestlambdaicl
	}
	parameters[["sigma"]] <- bestsig
	parametersicl[["mean"]] <- bestmeanicl
	parametersicl[["sigma"]] <- bestsigicl
	parameters[["df"]] <- dff1
	parametersicl[["df"]] <- dff1icl
	store[["fuzzy"]] <- bestz
	store[["parameters"]] <- parameters
	store[["allbic"]] <- bic[,Gs]
#	if(noconv){warning("At least one model did not run to convergence (most likely due to a non-invertible covariance matrix)")}
	icllist[["allicl"]] <- icl[,Gs]
	store[["bic"]] <- max(bic)
	icllist[["parameters"]] <- parametersicl
	icllist[["icl"]] <- max(icl)
	icllist[["fuzzy"]] <- bestzicl
	store[["bestmodel"]] <- paste("The best model (BIC of ",round(max(bic),2),") is ",bestmod," with G=",bestg,sep="")
	store[["classification"]] <- bestzmap
	icllist[["bestmodel"]] <- paste("The best model (ICL of ",round(max(icl),2),") is ",bestmodicl," with G=",bestgicl,sep="")
	icllist[["classification"]] <- bestzmapicl
	store[["G"]] <- bestg
	icllist[["G"]] <- bestgicl
	if(length(known)>0){
		store[["tab"]] <- tab
		icllist[["tab"]] <- tabicl
	}
	store[["x"]] <- x
	if(clas>0){
		store[["index"]] <- testindex
	}
	store[["logl"]] <- logls[which(bic==max(bic), arr.ind=TRUE)[1],which(bic==max(bic), arr.ind=TRUE)[2]]
	icllist[["logl"]] <- logls[which(icl==max(icl), arr.ind=TRUE)[1],which(icl==max(icl), arr.ind=TRUE)[2]]
	store[["iclresults"]] <- icllist
	store
}

